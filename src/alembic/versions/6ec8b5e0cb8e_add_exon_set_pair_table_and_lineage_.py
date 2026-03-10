"""add exon_set_pair table and lineage-aware views

Revision ID: 6ec8b5e0cb8e
Revises: 77076df4224c
Create Date: 2026-03-10 17:05:00.000000

"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = "6ec8b5e0cb8e"
down_revision: Union[str, None] = "77076df4224c"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.create_table(
        "exon_set_pair",
        sa.Column("alt_exon_set_id", sa.Integer(), nullable=False),
        sa.Column("tx_exon_set_id", sa.Integer(), nullable=False),
        sa.Column("added", sa.DateTime(), nullable=False),
        sa.ForeignKeyConstraint(
            ["alt_exon_set_id"],
            ["uta.exon_set.exon_set_id"],
            onupdate="CASCADE",
            ondelete="CASCADE",
        ),
        sa.ForeignKeyConstraint(
            ["tx_exon_set_id"],
            ["uta.exon_set.exon_set_id"],
            onupdate="CASCADE",
            ondelete="CASCADE",
        ),
        sa.PrimaryKeyConstraint("alt_exon_set_id"),
        schema="uta",
    )
    op.create_index(
        "ix_uta_exon_set_pair_tx_exon_set_id",
        "exon_set_pair",
        ["tx_exon_set_id"],
        unique=False,
        schema="uta",
    )

    op.execute("DROP MATERIALIZED VIEW IF EXISTS tx_exon_aln_mv CASCADE;")
    op.execute("DROP VIEW IF EXISTS tx_exon_aln_v CASCADE;")
    op.execute("DROP VIEW IF EXISTS tx_alt_exon_pairs_v CASCADE;")

    op.execute(
        """
        CREATE VIEW tx_alt_exon_pairs_v AS
        WITH effective_pairs AS (
            SELECT AES.exon_set_id AS aes_exon_set_id,
                   COALESCE(ESP.tx_exon_set_id, TES.exon_set_id) AS tes_exon_set_id,
                   CASE
                       WHEN ESP.tx_exon_set_id IS NOT NULL
                        AND AES.alt_aln_method LIKE '%/%'
                        AND NOT EXISTS (
                            SELECT 1
                            FROM exon_set plain
                            WHERE plain.tx_ac = AES.tx_ac
                              AND plain.alt_ac = AES.alt_ac
                              AND plain.alt_aln_method = regexp_replace(AES.alt_aln_method, '/.*$', '')
                        )
                       THEN regexp_replace(AES.alt_aln_method, '/.*$', '')
                       ELSE AES.alt_aln_method
                   END AS effective_alt_aln_method
            FROM exon_set AES
            JOIN exon_set TES
              ON TES.tx_ac = AES.tx_ac
             AND TES.alt_ac = AES.tx_ac
             AND TES.alt_aln_method = 'transcript'
            LEFT JOIN exon_set_pair ESP ON ESP.alt_exon_set_id = AES.exon_set_id
            WHERE AES.alt_aln_method !~ '^transcript(/|$)'
        )
        SELECT g.symbol, g.symbol as hgnc, g.gene_id, TES.exon_SET_id AS tes_exon_SET_id,
               AES.exon_SET_id AS aes_exon_SET_id, TES.tx_ac AS tx_ac, AES.alt_ac AS alt_ac,
               AES.alt_strand, EP.effective_alt_aln_method AS alt_aln_method, TEX.ORD AS ORD, TEX.exon_id AS tx_exon_id,
               AEX.exon_id AS alt_exon_id, TEX.start_i AS tx_start_i, TEX.END_i AS tx_END_i,
               AEX.start_i AS alt_start_i, AEX.END_i AS alt_END_i, EA.exon_aln_id, EA.cigar
        FROM effective_pairs EP
        JOIN exon_SET TES ON TES.exon_set_id = EP.tes_exon_set_id
        JOIN exon_SET AES ON AES.exon_set_id = EP.aes_exon_set_id
        JOIN transcript t ON TES.tx_ac = t.ac
        JOIN gene g ON t.gene_id = g.gene_id
        JOIN exon TEX ON TES.exon_SET_id = TEX.exon_SET_id
        JOIN exon AEX ON AES.exon_SET_id = AEX.exon_SET_id AND TEX.ORD = AEX.ORD
        LEFT JOIN exon_aln EA ON EA.tx_exon_id = TEX.exon_id AND EA.alt_exon_id = AEX.exon_id;
        """
    )

    op.execute(
        """
        CREATE VIEW tx_exon_aln_v AS
        WITH effective_pairs AS (
            SELECT AES.exon_set_id AS alt_exon_set_id,
                   AES.tx_ac,
                   COALESCE(ESP.tx_exon_set_id, TES.exon_set_id) AS tx_exon_set_id,
                   CASE
                       WHEN ESP.tx_exon_set_id IS NOT NULL
                        AND AES.alt_aln_method LIKE '%/%'
                        AND NOT EXISTS (
                            SELECT 1
                            FROM exon_set plain
                            WHERE plain.tx_ac = AES.tx_ac
                              AND plain.alt_ac = AES.alt_ac
                              AND plain.alt_aln_method = regexp_replace(AES.alt_aln_method, '/.*$', '')
                        )
                       THEN regexp_replace(AES.alt_aln_method, '/.*$', '')
                       ELSE AES.alt_aln_method
                   END AS effective_alt_aln_method
            FROM exon_set AES
            JOIN exon_set TES
              ON TES.tx_ac = AES.tx_ac
             AND TES.alt_ac = AES.tx_ac
             AND TES.alt_aln_method = 'transcript'
            LEFT JOIN exon_set_pair ESP ON ESP.alt_exon_set_id = AES.exon_set_id
            WHERE AES.alt_aln_method !~ '^transcript(/|$)'
        )
        SELECT g.symbol as hgnc, T.ac as tx_ac, AES.alt_ac, EP.effective_alt_aln_method AS alt_aln_method, AES.alt_strand,
               TE.ord, TE.start_i as tx_start_i, TE.end_i as tx_end_i,
               AE.start_i as alt_start_i, AE.end_i as alt_end_i,
               EA.cigar, EA.tx_aseq, EA.alt_aseq,
               TES.exon_set_id AS tx_exon_set_id, AES.exon_set_id as alt_exon_set_id,
               TE.exon_id as tx_exon_id, AE.exon_id as alt_exon_id,
               EA.exon_aln_id
        FROM effective_pairs EP
        JOIN transcript T ON EP.tx_ac = T.ac
        JOIN gene g ON T.gene_id = g.gene_id
        JOIN exon_set TES ON TES.exon_set_id = EP.tx_exon_set_id
        JOIN exon_set AES ON AES.exon_set_id = EP.alt_exon_set_id
        JOIN exon TE ON TES.exon_set_id = TE.exon_set_id
        JOIN exon AE ON AES.exon_set_id = AE.exon_set_id AND TE.ord = AE.ord
        LEFT JOIN exon_aln EA ON TE.exon_id = EA.tx_exon_id AND AE.exon_id = EA.alt_exon_id;
        """
    )

    op.execute(
        """
        CREATE MATERIALIZED VIEW tx_exon_aln_mv AS SELECT * FROM tx_exon_aln_v WITH NO DATA;
        GRANT SELECT ON tx_exon_aln_mv TO public;
        """
    )


def downgrade() -> None:
    op.execute("DROP MATERIALIZED VIEW IF EXISTS tx_exon_aln_mv CASCADE;")
    op.execute("DROP VIEW IF EXISTS tx_exon_aln_v CASCADE;")
    op.execute("DROP VIEW IF EXISTS tx_alt_exon_pairs_v CASCADE;")

    op.execute(
        """
        CREATE VIEW tx_alt_exon_pairs_v AS
            SELECT g.symbol, g.symbol as hgnc, g.gene_id, TES.exon_SET_id AS tes_exon_SET_id,
                   AES.exon_SET_id AS aes_exon_SET_id, TES.tx_ac AS tx_ac, AES.alt_ac AS alt_ac,
                   AES.alt_strand, AES.alt_aln_method, TEX.ORD AS ORD, TEX.exon_id AS tx_exon_id,
                   AEX.exon_id AS alt_exon_id, TEX.start_i AS tx_start_i, TEX.END_i AS tx_END_i,
                   AEX.start_i AS alt_start_i, AEX.END_i AS alt_END_i, EA.exon_aln_id, EA.cigar
            FROM exon_SET tes
            JOIN transcript t ON tes.tx_ac=t.ac
            JOIN gene g ON t.gene_id=g.gene_id
            JOIN exon_set aes ON tes.tx_ac=aes.tx_ac AND tes.alt_aln_method='transcript' AND aes.alt_aln_method !~ 'transcript'
            JOIN exon tex ON tes.exon_SET_id=tex.exon_SET_id
            JOIN exon aex ON aes.exon_SET_id=aex.exon_SET_id AND tex.ORD=aex.ORD
            LEFT JOIN exon_aln ea ON ea.tx_exon_id=tex.exon_id AND ea.alt_exon_id=AEX.exon_id;
        """
    )

    op.execute(
        """
        CREATE VIEW tx_exon_aln_v AS
            SELECT g.symbol as hgnc, T.ac as tx_ac, AES.alt_ac, AES.alt_aln_method, AES.alt_strand,
                   TE.ord, TE.start_i as tx_start_i, TE.end_i as tx_end_i,
                   AE.start_i as alt_start_i, AE.end_i as alt_end_i,
                   EA.cigar, EA.tx_aseq, EA.alt_aseq,
                   TES.exon_set_id AS tx_exon_set_id, AES.exon_set_id as alt_exon_set_id,
                   TE.exon_id as tx_exon_id, AE.exon_id as alt_exon_id,
                   EA.exon_aln_id
            FROM transcript T
            JOIN gene g ON T.gene_id = g.gene_id
            JOIN exon_set TES ON T.ac=TES.tx_ac AND TES.alt_aln_method='transcript'
            JOIN exon_set AES on T.ac=AES.tx_ac and AES.alt_aln_method !~ 'transcript'
            JOIN exon TE ON TES.exon_set_id=TE.exon_set_id
            JOIN exon AE ON AES.exon_set_id=AE.exon_set_id AND TE.ord=AE.ord
            LEFT JOIN exon_aln EA ON TE.exon_id=EA.tx_exon_id AND AE.exon_id=EA.alt_exon_id;
        """
    )

    op.execute(
        """
        CREATE MATERIALIZED VIEW tx_exon_aln_mv AS SELECT * FROM tx_exon_aln_v WITH NO DATA;
        GRANT SELECT ON tx_exon_aln_mv TO public;
        """
    )

    op.drop_index("ix_uta_exon_set_pair_tx_exon_set_id", table_name="exon_set_pair", schema="uta")
    op.drop_table("exon_set_pair", schema="uta")
