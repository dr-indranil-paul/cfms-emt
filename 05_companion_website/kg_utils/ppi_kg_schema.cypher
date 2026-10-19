// ============================================================
// schema_ppi.cypher — PPI-augmented Reactome knowledge graph
//
// Source: Reactome gk-central-06172025 (version 91)
//         + Emili Lab PPI data (fig2e_nodes.csv / fig2e_edges.csv)
//
// NOTE: dbId is intentionally NOT used as a primary unique constraint.
//   dbId is Reactome's internal integer and is not stable across versions.
//   Stable keys: stId (entities/events), identifier (reference nodes),
//   taxId (Species). dbId is indexed only as a lookup fallback for
//   metadata-only nodes (InstanceEdit, Person, Publication, etc.).
//
// Apply BEFORE importing ppi_kg_scaffold.cypher so that MERGE statements
// can use constraint-backed indexes rather than full scans.
//
// Run with:
//   docker exec <container> cypher-shell -u neo4j -p password < kg/schema_ppi.cypher
// ============================================================

// --- Stable keys for core Reactome node types ---

CREATE CONSTRAINT pathway_stid IF NOT EXISTS
    FOR (n:Pathway) REQUIRE n.stId IS UNIQUE;

CREATE CONSTRAINT rgp_identifier IF NOT EXISTS
    FOR (n:ReferenceGeneProduct) REQUIRE n.identifier IS UNIQUE;

CREATE CONSTRAINT go_bp_identifier IF NOT EXISTS
    FOR (n:GO_BiologicalProcess) REQUIRE n.identifier IS UNIQUE;

CREATE CONSTRAINT compartment_identifier IF NOT EXISTS
    FOR (n:Compartment) REQUIRE n.identifier IS UNIQUE;

CREATE CONSTRAINT species_taxid IF NOT EXISTS
    FOR (n:Species) REQUIRE n.taxId IS UNIQUE;

CREATE CONSTRAINT ref_therapeutic_identifier IF NOT EXISTS
    FOR (n:ReferenceTherapeutic) REQUIRE n.identifier IS UNIQUE;

// stId index on DatabaseObject for fast MERGE across events/entities
// (Pathway above is a unique constraint; this covers the broader set)
CREATE INDEX databaseobject_stid IF NOT EXISTS
    FOR (n:DatabaseObject) ON (n.stId);

// dbId index for metadata-only nodes that have no stId/identifier/taxId
// (InstanceEdit, Person, LiteratureReference, Summation, etc.)
CREATE INDEX databaseobject_dbid IF NOT EXISTS
    FOR (n:DatabaseObject) ON (n.dbId);

// variantIdentifier is the stable MERGE key for ReferenceIsoform nodes
// (isoforms share identifier with their canonical RGP — not usable as MERGE key)
CREATE INDEX reference_isoform_variant_id IF NOT EXISTS
    FOR (n:ReferenceIsoform) ON (n.variantIdentifier);

// --- PPI-added node types ---

// Assembly: one node per protein complex (C0001–C0506)
CREATE CONSTRAINT assembly_name IF NOT EXISTS
    FOR (n:Assembly) REQUIRE n.name IS UNIQUE;

// InteractionEvent: Reactome's native InteractionEvent nodes use stId (covered
// by databaseobject_stid above). PPI-added InteractionEvent nodes have no stId;
// they use name (format: "GENE1_GENE2_interaction"). No namespace collision.
CREATE INDEX interaction_event_name IF NOT EXISTS
    FOR (n:InteractionEvent) ON (n.name);

// ewas_category index speeds up MATCH on PPI EWAS nodes
// (PPI EWAS are merged on {name, ewas_category: 'ppi_ewas'})
CREATE INDEX ewas_ppi_category IF NOT EXISTS
    FOR (n:EntityWithAccessionedSequence) ON (n.ewas_category);

// --- Common query properties ---

CREATE INDEX rgp_display_name IF NOT EXISTS
    FOR (n:ReferenceGeneProduct) ON (n.displayName);

CREATE INDEX pathway_display_name IF NOT EXISTS
    FOR (n:Pathway) ON (n.displayName);

CREATE INDEX go_bp_display_name IF NOT EXISTS
    FOR (n:GO_BiologicalProcess) ON (n.displayName);
