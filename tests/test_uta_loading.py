import configparser
import gzip
import os
import signal
import tempfile
import unittest
from unittest.mock import Mock, patch

import sqlalchemy as sa
import testing.postgresql

import uta
import uta.loading as ul
import uta.models as usam


class TestUtaLoading(unittest.TestCase):

    def setUp(self):
        # setup test database
        self.db = testing.postgresql.Postgresql()
        self.session = uta.connect(self.db.url())

        admin_role = 'uta_admin'
        self.session.execute(sa.text(f'create user {admin_role}'))
        self.session.execute(sa.text(f'grant all privileges on database test to {admin_role}'))

        self.cf = configparser.ConfigParser()
        self.cf.add_section('uta')
        self.cf.set('uta', 'admin_role', 'uta_admin')

        ul.create_schema(self.session, {}, self.cf)
        ul.grant_permissions(self.session, {}, self.cf)

    def tearDown(self):
        self.session.close()
        self.db.stop(_signal=signal.SIGKILL)
        self.db.cleanup()

    def _create_exon_set_exons_view(self):
        self.session.execute(sa.text(
            """
            create or replace view uta.exon_set_exons_v as
            select ES.*,EL.n_exons,EL.se_i,EL.starts_i,EL.ends_i,EL.lengths
            from uta.exon_set ES
            join (
                select
                    iES.exon_set_id,
                    count(*) as n_exons,
                    array_to_string(array_agg(format('%s,%s', iE.start_i, iE.end_i) order by iE.ord), ';') as se_i,
                    array_agg(iE.start_i order by iE.ord) as starts_i,
                    array_agg(iE.end_i order by iE.ord) as ends_i,
                    array_agg((iE.end_i-iE.start_i) order by iE.ord) as lengths
                from uta.exon_set iES
                join uta.exon iE on iES.exon_set_id=iE.exon_set_id
                group by iES.exon_set_id
            ) EL on ES.exon_set_id = EL.exon_set_id
            """
        ))
        self.session.commit()

    def _create_lineage_aware_alignment_views(self):
        self.session.execute(sa.text(
            """
            create or replace view uta.tx_alt_exon_pairs_v as
            with effective_pairs as (
                select AES.exon_set_id as aes_exon_set_id,
                       coalesce(ESP.tx_exon_set_id, TES.exon_set_id) as tes_exon_set_id,
                       case
                           when ESP.tx_exon_set_id is not null
                            and AES.alt_aln_method like '%/%'
                            and not exists (
                                select 1
                                from uta.exon_set plain
                                where plain.tx_ac = AES.tx_ac
                                  and plain.alt_ac = AES.alt_ac
                                  and plain.alt_aln_method = regexp_replace(AES.alt_aln_method, '/.*$', '')
                            )
                           then regexp_replace(AES.alt_aln_method, '/.*$', '')
                           else AES.alt_aln_method
                       end as effective_alt_aln_method
                from uta.exon_set AES
                join uta.exon_set TES
                  on TES.tx_ac = AES.tx_ac
                 and TES.alt_ac = AES.tx_ac
                 and TES.alt_aln_method = 'transcript'
                left join uta.exon_set_pair ESP on ESP.alt_exon_set_id = AES.exon_set_id
                where AES.alt_aln_method !~ '^transcript(/|$)'
            )
            select g.symbol, g.symbol as hgnc, g.gene_id, TES.exon_set_id as tes_exon_set_id,
                   AES.exon_set_id as aes_exon_set_id, TES.tx_ac as tx_ac, AES.alt_ac as alt_ac,
                   AES.alt_strand, EP.effective_alt_aln_method as alt_aln_method, TEX.ord as ord, TEX.exon_id as tx_exon_id,
                   AEX.exon_id as alt_exon_id, TEX.start_i as tx_start_i, TEX.end_i as tx_end_i,
                   AEX.start_i as alt_start_i, AEX.end_i as alt_end_i, EA.exon_aln_id, EA.cigar
            from effective_pairs EP
            join uta.exon_set TES on TES.exon_set_id = EP.tes_exon_set_id
            join uta.exon_set AES on AES.exon_set_id = EP.aes_exon_set_id
            join uta.transcript T on TES.tx_ac = T.ac
            join uta.gene G on T.gene_id = G.gene_id
            join uta.exon TEX on TES.exon_set_id = TEX.exon_set_id
            join uta.exon AEX on AES.exon_set_id = AEX.exon_set_id and TEX.ord = AEX.ord
            left join uta.exon_aln EA on EA.tx_exon_id = TEX.exon_id and EA.alt_exon_id = AEX.exon_id
            """
        ))
        self.session.execute(sa.text(
            """
            create or replace view uta.tx_exon_aln_v as
            with effective_pairs as (
                select AES.exon_set_id as alt_exon_set_id,
                       AES.tx_ac,
                       coalesce(ESP.tx_exon_set_id, TES.exon_set_id) as tx_exon_set_id,
                       case
                           when ESP.tx_exon_set_id is not null
                            and AES.alt_aln_method like '%/%'
                            and not exists (
                                select 1
                                from uta.exon_set plain
                                where plain.tx_ac = AES.tx_ac
                                  and plain.alt_ac = AES.alt_ac
                                  and plain.alt_aln_method = regexp_replace(AES.alt_aln_method, '/.*$', '')
                            )
                           then regexp_replace(AES.alt_aln_method, '/.*$', '')
                           else AES.alt_aln_method
                       end as effective_alt_aln_method
                from uta.exon_set AES
                join uta.exon_set TES
                  on TES.tx_ac = AES.tx_ac
                 and TES.alt_ac = AES.tx_ac
                 and TES.alt_aln_method = 'transcript'
                left join uta.exon_set_pair ESP on ESP.alt_exon_set_id = AES.exon_set_id
                where AES.alt_aln_method !~ '^transcript(/|$)'
            )
            select G.symbol as hgnc, T.ac as tx_ac, AES.alt_ac,
                   EP.effective_alt_aln_method as alt_aln_method, AES.alt_strand,
                   TE.ord, TE.start_i as tx_start_i, TE.end_i as tx_end_i,
                   AE.start_i as alt_start_i, AE.end_i as alt_end_i,
                   EA.cigar, EA.tx_aseq, EA.alt_aseq,
                   TES.exon_set_id as tx_exon_set_id, AES.exon_set_id as alt_exon_set_id,
                   TE.exon_id as tx_exon_id, AE.exon_id as alt_exon_id,
                   EA.exon_aln_id
            from effective_pairs EP
            join uta.transcript T on EP.tx_ac = T.ac
            join uta.gene G on T.gene_id = G.gene_id
            join uta.exon_set TES on TES.exon_set_id = EP.tx_exon_set_id
            join uta.exon_set AES on AES.exon_set_id = EP.alt_exon_set_id
            join uta.exon TE on TES.exon_set_id = TE.exon_set_id
            join uta.exon AE on AES.exon_set_id = AE.exon_set_id and TE.ord = AE.ord
            left join uta.exon_aln EA on TE.exon_id = EA.tx_exon_id and AE.exon_id = EA.alt_exon_id
            """
        ))
        self.session.commit()

    def test_meta_data(self):
        """
        Metadata should exist, then updated when update_meta_data is called.
        """
        # the schema_version should match existing values in UTA models
        expected_schema_version = usam.schema_version
        md_schema_version = self.session.query(usam.Meta).filter(usam.Meta.key == 'schema_version').one()
        self.assertEqual(md_schema_version.value, expected_schema_version)

        new_schema_version = '9.9'
        with patch('uta.models.schema_version', new_schema_version):
            ul.update_meta_data(self.session, {}, self.cf)

        md_schema_version = self.session.query(usam.Meta).filter(usam.Meta.key == 'schema_version').one()
        self.assertEqual(md_schema_version.value, new_schema_version)

        md_updated_at = self.session.query(usam.Meta).filter(usam.Meta.key == 'updated on').one_or_none()
        self.assertIsNotNone(md_updated_at)

    def test_load_assoc_ac(self):
        """
        Loading file tests/data/assocacs.gz should create associated_accessions records in the database.
        Row will be created in associated_accessions even when transcript or origin does not exist in database.
        This is only the case until tx_ac and origin are converted to foreign keys.
        """

        # insert origins referenced in data file
        o1 = usam.Origin(
            name='NCBI',
            url='http://bogus.com/ncbi',
            url_ac_fmt='http://bogus.com/ncbi/{ac}',
        )
        self.session.add(o1)

        # insert genes required for transcripts
        g1 = usam.Gene(
            gene_id='49',
            hgnc='ACR',
            symbol='ACR',
            maploc='22q13.33',
            descr='acrosin',
            summary='acrosin',
            aliases='SPGF87',
            type='protein-coding',
            xrefs='MIM:102480,HGNC:HGNC:126,Ensembl:ENSG00000100312,AllianceGenome:HGNC:126',
        )
        g2 = usam.Gene(
            gene_id=50,
            hgnc='ACO2',
            symbol='ACO2',
            maploc='22q13.2',
            descr='aconitase 2',
            summary='aconitase 2',
            aliases='ACONM,HEL-S-284,ICRD,OCA8,OPA9',
            type='protein-coding',
            xrefs='MIM:100850,HGNC:HGNC:118,Ensembl:ENSG00000100412,AllianceGenome:HGNC:118',
        )
        self.session.add(g1)
        self.session.add(g2)

        # insert transcripts referenced in data file
        t1 = usam.Transcript(
            ac='NM_001097.3',
            origin=o1,
            gene_id=g1.gene_id,
            cds_start_i=0,
            cds_end_i=1,
            cds_md5='a',
        )
        t2 = usam.Transcript(
            ac='NM_001098.3',
            origin=o1,
            gene_id=g2.gene_id,
            cds_start_i=2,
            cds_end_i=3,
            cds_md5='b',
        )
        self.session.add(t1)
        self.session.add(t2)

        # pre-add one of the associated_acessions from the test data file
        # to demonstrate get-or-insert behavior
        p = usam.AssociatedAccessions(
            tx_ac='NM_001097.3',
            pro_ac='NP_001088.2',
            origin='NCBI',
        )
        self.session.add(p)

        self.session.commit()

        cf = configparser.ConfigParser()
        cf.add_section('uta')
        cf.set('uta', 'admin_role', 'uta_admin')

        ul.load_assoc_ac(self.session, {'FILE': 'tests/data/assocacs.gz'}, cf)

        # associated_accessions table should contain one record per line in file
        aa = self.session.query(usam.AssociatedAccessions).order_by(usam.AssociatedAccessions.tx_ac).all()
        aa_list = [{'tx_ac': aa.tx_ac, 'pro_ac': aa.pro_ac, 'origin_name': aa.origin} for aa in aa]
        expected_aa_list = [
            {
                'tx_ac': 'DummyTx',
                'pro_ac': 'DummyProtein',
                'origin_name': 'DummyOrigin',
            },
            {
                'tx_ac': 'NM_001097.3',
                'pro_ac': 'NP_001088.2',
                'origin_name': 'NCBI',
            },
            {
                'tx_ac': 'NM_001098.3',
                'pro_ac': 'NP_001089.1',
                'origin_name': 'NCBI',
            },
        ]
        self.assertEqual(aa_list, expected_aa_list)

    def test_load_txinfo(self):
        """
        Loading file tests/data/txinfo.gz should create transcript, exon_set, exon, and translation_exception records in the database.
        """

        # insert origins referenced in data file
        o1 = usam.Origin(
            name='NCBI',
            url='http://bogus.com/ncbi',
            url_ac_fmt='http://bogus.com/ncbi/{ac}',
        )
        self.session.add(o1)

        # insert genes required for transcripts
        g1 = usam.Gene(
            gene_id='140606',
            hgnc='SELENOM',
            symbol='SELENOM',
            maploc='22q12.2',
            descr='selenoprotein M',
            summary='selenoprotein M',
            aliases='SELM,SEPM',
            type='protein-coding',
            xrefs='MIM:610918,HGNC:HGNC:30397,Ensembl:ENSG00000198832,AllianceGenome:HGNC:30397',
        )
        g2 = usam.Gene(
            gene_id='4514',
            hgnc='MT-CO3',
            symbol='MT-CO3',
            maploc=None,
            descr='mitochondrially encoded cytochrome c oxidase III',
            summary='mitochondrially encoded cytochrome c oxidase III',
            aliases='COIII,MTCO3',
            type='protein-coding',
            xrefs='GeneID:4514,HGNC:HGNC:7422,MIM:516050',
        )
        self.session.add(g1)
        self.session.add(g2)
        self.session.commit()

        cf = configparser.ConfigParser()
        cf.add_section('uta')
        cf.set('uta', 'admin_role', 'uta_admin')

        with patch('uta.loading._get_seqfetcher', return_value=Mock(fetch=Mock(return_value='FAKESEQUENCE'))):
            ul.load_txinfo(self.session, {'FILE': 'tests/data/txinfo.gz'}, cf)

        transcript = self.session.query(usam.Transcript).filter(usam.Transcript.ac == 'NM_080430.4').one()
        self.assertEqual(
            {
                'ac': transcript.ac,
                'gene_id': transcript.gene_id,
                'cds_start_i': transcript.cds_start_i,
                'cds_end_i': transcript.cds_end_i,
                'codon_table': transcript.codon_table,
            },
            {
                'ac': 'NM_080430.4',
                'gene_id': '140606',
                'cds_start_i': 63,
                'cds_end_i': 501,
                'codon_table': '1',
            },
        )

        transcript = self.session.query(usam.Transcript).filter(usam.Transcript.ac == 'NC_012920.1_09206_09990').one()
        self.assertEqual(
            {
                'ac': transcript.ac,
                'gene_id': transcript.gene_id,
                'cds_start_i': transcript.cds_start_i,
                'cds_end_i': transcript.cds_end_i,
                'codon_table': transcript.codon_table,
            },
            {
                'ac': 'NC_012920.1_09206_09990',
                'gene_id': '4514',
                'cds_start_i': 0,
                'cds_end_i': 784,
                'codon_table': '2',
            },
        )

        exon_set = self.session.query(usam.ExonSet).filter(usam.ExonSet.tx_ac == 'NM_080430.4').one()
        exons = self.session.query(usam.Exon).filter(usam.Exon.exon_set_id == exon_set.exon_set_id).all()
        self.assertEqual(len(exons), 5)

        translation_exception = self.session.query(usam.TranslationException).filter(usam.TranslationException.tx_ac == 'NM_080430.4').one()
        self.assertEqual(
            {
                'tx_ac': translation_exception.tx_ac,
                'start_position': translation_exception.start_position,
                'end_position': translation_exception.end_position,
                'amino_acid': translation_exception.amino_acid,
            },
            {
                'tx_ac': 'NM_080430.4',
                'start_position': 204,
                'end_position': 207,
                'amino_acid': 'Sec',
            },
        )

    def test_load_txinfo_archives_plain_alt_exon_sets_when_transcript_exons_change(self):
        o1 = usam.Origin(
            name='NCBI',
            url='http://bogus.com/ncbi',
            url_ac_fmt='http://bogus.com/ncbi/{ac}',
        )
        self.session.add(o1)

        g1 = usam.Gene(
            gene_id='1',
            hgnc='GENE1',
            symbol='GENE1',
            maploc='1q1.1',
            descr='test gene',
            summary='test gene',
            aliases='GENE1',
            type='protein-coding',
            xrefs='GeneID:1',
        )
        self.session.add(g1)
        self.session.flush()

        tx = usam.Transcript(
            ac='NM_TEST.1',
            origin_id=o1.origin_id,
            gene_id=g1.gene_id,
            cds_start_i=0,
            cds_end_i=20,
            cds_md5='a',
            codon_table='1',
        )
        self.session.add(tx)
        self.session.flush()

        tx_es = usam.ExonSet(
            tx_ac=tx.ac,
            alt_ac=tx.ac,
            alt_strand=1,
            alt_aln_method='transcript',
        )
        splign_es = usam.ExonSet(
            tx_ac=tx.ac,
            alt_ac='NC_000001.11',
            alt_strand=1,
            alt_aln_method='splign',
        )
        blat_es = usam.ExonSet(
            tx_ac=tx.ac,
            alt_ac='AC_000001.1',
            alt_strand=1,
            alt_aln_method='blat',
        )
        self.session.add_all([tx_es, splign_es, blat_es])
        self.session.flush()

        for exon_set, exons in [
            (tx_es, [(0, 10), (10, 20)]),
            (splign_es, [(100, 110), (110, 120)]),
            (blat_es, [(200, 210), (210, 220)]),
        ]:
            for ord_, (start_i, end_i) in enumerate(exons):
                self.session.add(
                    usam.Exon(
                        exon_set_id=exon_set.exon_set_id,
                        start_i=start_i,
                        end_i=end_i,
                        ord=ord_,
                    )
                )
        self.session.commit()

        fd, txinfo_path = tempfile.mkstemp(suffix='.gz')
        os.close(fd)
        try:
            with gzip.open(txinfo_path, 'wt') as f:
                f.write("origin\tac\tgene_id\tgene_symbol\tcds_se_i\texons_se_i\tcodon_table\ttransl_except\n")
                f.write("NCBI\tNM_TEST.1\t1\tGENE1\t0,20\t0,10;10,20;20,30\t1\t\n")

            with patch('uta.loading._get_seqfetcher', return_value=Mock(fetch=Mock(return_value='FAKESEQUENCE'))):
                ul.load_txinfo(self.session, {'FILE': txinfo_path}, self.cf)
        finally:
            os.unlink(txinfo_path)

        plain_tx_es = self.session.query(usam.ExonSet).filter(
            usam.ExonSet.tx_ac == 'NM_TEST.1',
            usam.ExonSet.alt_ac == 'NM_TEST.1',
            usam.ExonSet.alt_aln_method == 'transcript',
        ).one()
        self.assertEqual(len(plain_tx_es.exons), 3)

        archived_tx_es = self.session.query(usam.ExonSet).filter(
            usam.ExonSet.tx_ac == 'NM_TEST.1',
            usam.ExonSet.alt_ac == 'NM_TEST.1',
            usam.ExonSet.alt_aln_method.like('transcript/%'),
        ).one()
        self.assertEqual(len(archived_tx_es.exons), 2)

        for alt_ac, method in [('NC_000001.11', 'splign'), ('AC_000001.1', 'blat')]:
            plain_alt_count = self.session.query(usam.ExonSet).filter(
                usam.ExonSet.tx_ac == 'NM_TEST.1',
                usam.ExonSet.alt_ac == alt_ac,
                usam.ExonSet.alt_aln_method == method,
            ).count()
            self.assertEqual(plain_alt_count, 0)

            archived_alt_es = self.session.query(usam.ExonSet).filter(
                usam.ExonSet.tx_ac == 'NM_TEST.1',
                usam.ExonSet.alt_ac == alt_ac,
                usam.ExonSet.alt_aln_method.like(method + '/%'),
            ).one()
            self.assertEqual(len(archived_alt_es.exons), 2)

    def test_archive_exon_set_record_deletes_plain_duplicate_when_hashed_row_exists(self):
        o1 = usam.Origin(
            name='NCBI',
            url='http://bogus.com/ncbi',
            url_ac_fmt='http://bogus.com/ncbi/{ac}',
        )
        self.session.add(o1)

        g1 = usam.Gene(
            gene_id='1',
            hgnc='GENE1',
            symbol='GENE1',
            maploc='1q1.1',
            descr='test gene',
            summary='test gene',
            aliases='GENE1',
            type='protein-coding',
            xrefs='GeneID:1',
        )
        self.session.add(g1)
        self.session.flush()

        tx = usam.Transcript(
            ac='NM_DUP.1',
            origin_id=o1.origin_id,
            gene_id='1',
            cds_start_i=0,
            cds_end_i=20,
            cds_md5='a',
        )
        self.session.add(tx)
        self.session.flush()

        plain_es = usam.ExonSet(
            tx_ac=tx.ac,
            alt_ac='NC_000001.11',
            alt_strand=1,
            alt_aln_method='splign',
        )
        archived_es = usam.ExonSet(
            tx_ac=tx.ac,
            alt_ac='NC_000001.11',
            alt_strand=1,
            alt_aln_method='splign/9c2a62cb',
        )
        self.session.add_all([plain_es, archived_es])
        self.session.flush()

        for exon_set in [plain_es, archived_es]:
            for ord_, (start_i, end_i) in enumerate([(100, 110), (110, 120), (120, 130)]):
                self.session.add(
                    usam.Exon(
                        exon_set_id=exon_set.exon_set_id,
                        start_i=start_i,
                        end_i=end_i,
                        ord=ord_,
                    )
                )
        self.session.commit()

        result = ul._archive_exon_set_record(self.session, plain_es)
        self.session.commit()

        self.assertEqual(result.exon_set_id, archived_es.exon_set_id)
        self.assertEqual(
            self.session.query(usam.ExonSet).filter(
                usam.ExonSet.tx_ac == 'NM_DUP.1',
                usam.ExonSet.alt_ac == 'NC_000001.11',
                usam.ExonSet.alt_aln_method == 'splign',
            ).count(),
            0,
        )
        surviving = self.session.query(usam.ExonSet).filter(
            usam.ExonSet.tx_ac == 'NM_DUP.1',
            usam.ExonSet.alt_ac == 'NC_000001.11',
            usam.ExonSet.alt_aln_method == 'splign/9c2a62cb',
        ).one()
        self.assertEqual(len(surviving.exons), 3)

    def test_repair_stale_exon_set_overlays_archives_mismatches_and_promotes_matches(self):
        self._create_exon_set_exons_view()

        o1 = usam.Origin(
            name='NCBI',
            url='http://bogus.com/ncbi',
            url_ac_fmt='http://bogus.com/ncbi/{ac}',
        )
        self.session.add(o1)

        g1 = usam.Gene(
            gene_id='1',
            hgnc='GENE1',
            symbol='GENE1',
            maploc='1q1.1',
            descr='test gene',
            summary='test gene',
            aliases='GENE1',
            type='protein-coding',
            xrefs='GeneID:1',
        )
        self.session.add(g1)
        self.session.flush()

        tx = usam.Transcript(
            ac='NM_REPAIR.1',
            origin_id=o1.origin_id,
            gene_id=g1.gene_id,
            cds_start_i=0,
            cds_end_i=30,
            cds_md5='a',
            codon_table='1',
        )
        self.session.add(tx)
        self.session.flush()

        exon_sets = [
            ('NM_REPAIR.1', 'NM_REPAIR.1', 'transcript', [(0, 10), (10, 20), (20, 30)]),
            ('NM_REPAIR.1', 'NM_REPAIR.1', 'transcript/oldhash', [(0, 15), (15, 30)]),
            ('NM_REPAIR.1', 'NC_GOOD', 'splign', [(100, 110), (110, 120), (120, 130)]),
            ('NM_REPAIR.1', 'NC_BAD', 'splign', [(200, 215), (215, 230)]),
            ('NM_REPAIR.1', 'NC_PROMOTE', 'splign/newhash', [(300, 310), (310, 320), (320, 330)]),
        ]

        for tx_ac, alt_ac, method, exons in exon_sets:
            es = usam.ExonSet(
                tx_ac=tx_ac,
                alt_ac=alt_ac,
                alt_strand=1,
                alt_aln_method=method,
            )
            self.session.add(es)
            self.session.flush()
            for ord_, (start_i, end_i) in enumerate(exons):
                self.session.add(
                    usam.Exon(
                        exon_set_id=es.exon_set_id,
                        start_i=start_i,
                        end_i=end_i,
                        ord=ord_,
                    )
                )
        self.session.commit()
        ul.refresh_matviews(self.session, {}, self.cf)

        mismatch_rows = ul._plain_alt_exon_count_mismatch_rows(self.session)
        self.assertEqual(len(mismatch_rows), 1)
        self.assertEqual(mismatch_rows[0]["alt_ac"], "NC_BAD")

        with self.assertRaises(RuntimeError):
            ul.validate_exon_set_consistency(self.session, {}, self.cf)

        ul.repair_stale_exon_set_overlays(self.session, {}, self.cf)

        ul.validate_exon_set_consistency(self.session, {}, self.cf)

        self.assertEqual(
            self.session.query(usam.ExonSet).filter(
                usam.ExonSet.tx_ac == 'NM_REPAIR.1',
                usam.ExonSet.alt_ac == 'NC_GOOD',
                usam.ExonSet.alt_aln_method == 'splign',
            ).count(),
            1,
        )
        self.assertEqual(
            self.session.query(usam.ExonSet).filter(
                usam.ExonSet.tx_ac == 'NM_REPAIR.1',
                usam.ExonSet.alt_ac == 'NC_BAD',
                usam.ExonSet.alt_aln_method == 'splign',
            ).count(),
            1,
        )
        self.assertEqual(
            self.session.query(usam.ExonSet).filter(
                usam.ExonSet.tx_ac == 'NM_REPAIR.1',
                usam.ExonSet.alt_ac == 'NC_BAD',
                usam.ExonSet.alt_aln_method.like('splign/%'),
            ).count(),
            0,
        )
        self.assertEqual(
            self.session.query(usam.ExonSet).filter(
                usam.ExonSet.tx_ac == 'NM_REPAIR.1',
                usam.ExonSet.alt_ac == 'NC_PROMOTE',
                usam.ExonSet.alt_aln_method == 'splign',
            ).count(),
            1,
        )
        self.assertEqual(
            self.session.query(usam.ExonSet).filter(
                usam.ExonSet.tx_ac == 'NM_REPAIR.1',
                usam.ExonSet.alt_ac == 'NC_PROMOTE',
                usam.ExonSet.alt_aln_method.like('splign/%'),
            ).count(),
            0,
        )
        self.assertEqual(
            self.session.query(usam.ExonSetPair).filter(
                usam.ExonSetPair.alt_exon_set.has(tx_ac='NM_REPAIR.1'),
                usam.ExonSetPair.alt_exon_set.has(alt_ac='NC_BAD'),
            ).count(),
            1,
        )

    def test_pair_historical_exon_sets_pairs_archived_alt_to_matching_transcript_lineage(self):
        self._create_exon_set_exons_view()

        o1 = usam.Origin(
            name='NCBI',
            url='http://bogus.com/ncbi',
            url_ac_fmt='http://bogus.com/ncbi/{ac}',
        )
        self.session.add(o1)

        g1 = usam.Gene(
            gene_id='1',
            hgnc='GENE1',
            symbol='GENE1',
            maploc='1q1.1',
            descr='test gene',
            summary='test gene',
            aliases='GENE1',
            type='protein-coding',
            xrefs='GeneID:1',
        )
        self.session.add(g1)
        self.session.flush()

        tx = usam.Transcript(
            ac='NM_PAIR.1',
            origin_id=o1.origin_id,
            gene_id=g1.gene_id,
            cds_start_i=0,
            cds_end_i=30,
            cds_md5='a',
            codon_table='1',
        )
        self.session.add(tx)
        self.session.flush()

        exon_sets = [
            ('NM_PAIR.1', 'NM_PAIR.1', 'transcript', [(0, 10), (10, 20), (20, 30)]),
            ('NM_PAIR.1', 'NM_PAIR.1', 'transcript/oldhash', [(0, 15), (15, 30)]),
            ('NM_PAIR.1', 'NC_MATCH', 'splign/oldhash', [(100, 115), (115, 130)]),
            ('NM_PAIR.1', 'NC_PLAIN', 'splign', [(200, 210), (210, 220), (220, 230)]),
        ]

        created = {}
        for tx_ac, alt_ac, method, exons in exon_sets:
            es = usam.ExonSet(
                tx_ac=tx_ac,
                alt_ac=alt_ac,
                alt_strand=1,
                alt_aln_method=method,
            )
            self.session.add(es)
            self.session.flush()
            created[(alt_ac, method)] = es
            for ord_, (start_i, end_i) in enumerate(exons):
                self.session.add(
                    usam.Exon(
                        exon_set_id=es.exon_set_id,
                        start_i=start_i,
                        end_i=end_i,
                        ord=ord_,
                    )
                )
        self.session.commit()

        ul.pair_historical_exon_sets(self.session, {}, self.cf)
        ul.validate_historical_exon_set_pairings(self.session, {}, self.cf)

        pair = self.session.query(usam.ExonSetPair).filter(
            usam.ExonSetPair.alt_exon_set_id == created[('NC_MATCH', 'splign/oldhash')].exon_set_id
        ).one()
        self.assertEqual(
            pair.tx_exon_set_id,
            created[('NM_PAIR.1', 'transcript/oldhash')].exon_set_id,
        )
        self.assertEqual(
            self.session.query(usam.ExonSet).filter(
                usam.ExonSet.exon_set_id == created[('NC_MATCH', 'splign/oldhash')].exon_set_id
            ).one().alt_aln_method,
            'splign',
        )

    def test_lineage_aware_views_use_paired_transcript_lineage(self):
        self._create_exon_set_exons_view()
        self._create_lineage_aware_alignment_views()

        o1 = usam.Origin(
            name='NCBI',
            url='http://bogus.com/ncbi',
            url_ac_fmt='http://bogus.com/ncbi/{ac}',
        )
        self.session.add(o1)

        g1 = usam.Gene(
            gene_id='1',
            hgnc='GENE1',
            symbol='GENE1',
            maploc='1q1.1',
            descr='test gene',
            summary='test gene',
            aliases='GENE1',
            type='protein-coding',
            xrefs='GeneID:1',
        )
        self.session.add(g1)
        self.session.flush()

        tx = usam.Transcript(
            ac='NM_VIEW.1',
            origin_id=o1.origin_id,
            gene_id=g1.gene_id,
            cds_start_i=0,
            cds_end_i=30,
            cds_md5='a',
            codon_table='1',
        )
        self.session.add(tx)
        self.session.flush()

        exon_sets = [
            ('NM_VIEW.1', 'NM_VIEW.1', 'transcript', [(0, 10), (10, 20), (20, 30)]),
            ('NM_VIEW.1', 'NM_VIEW.1', 'transcript/oldhash', [(0, 15), (15, 30)]),
            ('NM_VIEW.1', 'NC_OLD', 'splign/oldhash', [(100, 115), (115, 130)]),
            ('NM_VIEW.1', 'NC_NEW', 'splign', [(200, 210), (210, 220), (220, 230)]),
        ]

        created = {}
        for tx_ac, alt_ac, method, exons in exon_sets:
            es = usam.ExonSet(
                tx_ac=tx_ac,
                alt_ac=alt_ac,
                alt_strand=1,
                alt_aln_method=method,
            )
            self.session.add(es)
            self.session.flush()
            created[(alt_ac, method)] = es
            for ord_, (start_i, end_i) in enumerate(exons):
                self.session.add(
                    usam.Exon(
                        exon_set_id=es.exon_set_id,
                        start_i=start_i,
                        end_i=end_i,
                        ord=ord_,
                    )
                )
        self.session.commit()

        self.session.add(
            usam.ExonSetPair(
                alt_exon_set_id=created[('NC_OLD', 'splign/oldhash')].exon_set_id,
                tx_exon_set_id=created[('NM_VIEW.1', 'transcript/oldhash')].exon_set_id,
            )
        )
        self.session.commit()

        old_rows = list(self.session.execute(sa.text(
            """
            select distinct tx_exon_set_id, alt_exon_set_id, alt_aln_method
            from uta.tx_exon_aln_v
            where tx_ac = 'NM_VIEW.1' and alt_ac = 'NC_OLD'
            """
        )))
        self.assertEqual(len(old_rows), 1)
        self.assertEqual(old_rows[0].tx_exon_set_id, created[('NM_VIEW.1', 'transcript/oldhash')].exon_set_id)
        self.assertEqual(old_rows[0].alt_exon_set_id, created[('NC_OLD', 'splign/oldhash')].exon_set_id)
        self.assertEqual(old_rows[0].alt_aln_method, 'splign')

        new_rows = list(self.session.execute(sa.text(
            """
            select distinct tx_exon_set_id, alt_exon_set_id, alt_aln_method
            from uta.tx_exon_aln_v
            where tx_ac = 'NM_VIEW.1' and alt_ac = 'NC_NEW'
            """
        )))
        self.assertEqual(len(new_rows), 1)
        self.assertEqual(new_rows[0].tx_exon_set_id, created[('NM_VIEW.1', 'transcript')].exon_set_id)
        self.assertEqual(new_rows[0].alt_exon_set_id, created[('NC_NEW', 'splign')].exon_set_id)
        self.assertEqual(new_rows[0].alt_aln_method, 'splign')

    def test_load_exonset_with_exon_structure_mismatch(self):
        """
        Loading the test file tests/data/exonsets-mm-exons.gz should not raise an exception, exon alignments without
        a mismatch should load, those with a mismatch should be skipped and logged as such. The input file has
        alignments for 4 transcripts against NC_000001.11, but only 2 of them have the correct number of exons.
        We only expect the alignmets for NM_000911.4 and NM_001005277.1 to be loaded.
        """
        # setup
        # insert origins referenced in data file
        o1 = usam.Origin(
            name="NCBI",
            url="http://bogus.com/ncbi",
            url_ac_fmt="http://bogus.com/ncbi/{ac}",
        )
        self.session.add(o1)
        self.session.flush()

        for gene_data in [
            {
                "gene_id": "3352",
                "hgnc": "HTR1D",
                "symbol": "HTR1D",
                "type": "protein-coding",
            },
            {
                "gene_id": "4985",
                "hgnc": "OPRD1",
                "symbol": "OPRD1",
                "type": "protein-coding",
            },
            {
                "gene_id": "81399",
                "hgnc": "OR4F16",
                "symbol": "OR4F16",
                "type": "protein-coding",
            },
            {
                "gene_id": "79501",
                "hgnc": "OR4F5",
                "symbol": "OR4F5",
                "type": "protein-coding",
            },
        ]:
            gene = usam.Gene(**gene_data)
            self.session.add(gene)

        for tx_data in [
            {
                "ac": "NM_000864.5",
                "origin_id": o1.origin_id,
                "gene_id": "3352",
                "cds_start_i": 994,
                "cds_end_i": 2128,
                "cds_md5": "a",
            },
            {
                "ac": "NM_000911.4",
                "origin_id": o1.origin_id,
                "gene_id": "4985",
                "cds_start_i": 214,
                "cds_end_i": 1333,
                "cds_md5": "b",
            },
            {
                "ac": "NM_001005277.1",
                "origin_id": o1.origin_id,
                "gene_id": "81399",
                "cds_start_i": 0,
                "cds_end_i": 939,
                "cds_md5": "c",
            },
            {
                "ac": "NM_001005484.2",
                "origin_id": o1.origin_id,
                "gene_id": "79501",
                "cds_start_i": 60,
                "cds_end_i": 1041,
                "cds_md5": "d",
            },
        ]:
            tx = usam.Transcript(**tx_data)
            self.session.add(tx)
            es = usam.ExonSet(
                tx_ac=tx.ac,
                alt_ac=tx.ac,
                alt_strand=1,
                alt_aln_method="transcript",
            )
            self.session.add(es)
            self.session.flush()

        for exon_data in [
            ("NM_000864.5", 1, 0, 3319),  # exons for NM_000864.5 are 0,212;212,3319
            ("NM_000911.4", 1, 0, 441),
            ("NM_000911.4", 2, 441, 791),
            ("NM_000911.4", 3, 791, 9317),
            ("NM_001005277.1", 1, 0, 939),
            ("NM_001005484.2", 1, 0, 15),
            ("NM_001005484.2", 2, 15, 69),
            (
                "NM_001005484.2",
                3,
                69,
                1041,
            ),  # exons for NM_001005484.2 are 0,15;15,69;69,2618
            ("NM_001005484.2", 4, 1041, 2618),
        ]:
            es = (
                self.session.query(usam.ExonSet)
                .filter(
                    usam.ExonSet.tx_ac == exon_data[0], usam.ExonSet.alt_ac == exon_data[0]
                )
                .one()
            )
            exon = usam.Exon(
                exon_set_id=es.exon_set_id,
                start_i=exon_data[2],
                end_i=exon_data[3],
                ord=exon_data[1],
            )
            self.session.add(exon)
        self.session.commit()

        cf = configparser.ConfigParser()
        cf.add_section("uta")
        cf.set("uta", "admin_role", "uta_admin")

        # load data from test exonsets file.
        with patch(
            "uta.loading._get_seqfetcher",
            return_value=Mock(fetch=Mock(return_value="FAKESEQUENCE")),
        ), patch("uta.loading.logger") as mock_logger:
            ul.load_exonset(self.session, {"FILE": "tests/data/exonsets.mm-exons.gz"}, cf)

            assert mock_logger.warning.called_with(
                "Exon structure mismatch: 4 exons in transcript NM_001005484.2; 3 in alignment NC_000001.11"
            )
            assert mock_logger.warning.called_with(
                "Exon structure mismatch: 1 exons in transcript NM_000864.5; 2 in alignment NC_000001.11"
            )

        # check that the exons for NM_000864.5 and NM_001005484.2 were not loaded,
        # and NM_000911.4 and NM_001005277.1 were loaded
        for tx_ac, expected_exon_count in [("NM_000911.4", 3), ("NM_001005277.1", 1)]:
            exon_set = (
                self.session.query(usam.ExonSet)
                .filter(
                    usam.ExonSet.tx_ac == tx_ac,
                    usam.ExonSet.alt_ac == "NC_000001.11",
                    usam.ExonSet.alt_aln_method == "splign",
                )
                .one()
            )
            exons = (
                self.session.query(usam.Exon)
                .filter(usam.Exon.exon_set_id == exon_set.exon_set_id)
                .all()
            )
            self.assertEqual(len(exons), expected_exon_count)

        for tx_ac in ["NM_000864.5", "NM_001005484.2"]:
            with self.assertRaises(sa.orm.exc.NoResultFound):
                self.session.query(usam.ExonSet).filter(
                    usam.ExonSet.tx_ac == tx_ac,
                    usam.ExonSet.alt_ac == "NC_000001.11",
                    usam.ExonSet.alt_aln_method == "splign",
                ).one()


class TestUtaLoadingFunctions(unittest.TestCase):
    def test__create_translation_exceptions(self):
        transl_except_list = ['(pos:333..335,aa:Sec)', '(pos:1017,aa:TERM)']
        translation_exceptions = ul._create_translation_exceptions(transcript='dummy_tx', transl_except_list=transl_except_list)
        self.assertEqual(translation_exceptions, [
            {
                'tx_ac': 'dummy_tx',
                'start_position': 332,
                'end_position': 335,
                'amino_acid': 'Sec',
            },
            {
                'tx_ac': 'dummy_tx',
                'start_position': 1016,
                'end_position': 1017,
                'amino_acid': 'TERM',
            },
        ])
