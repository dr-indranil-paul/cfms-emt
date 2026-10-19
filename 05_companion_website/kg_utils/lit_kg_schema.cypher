// ============================================================
// schema.cypher — ChatEMT literature knowledge graph
//
// Source: Reactome DB gk-central-06172025 (version 91)
//
// Apply BEFORE importing scaffold.cypher so that MERGE statements
// can use constraint-backed indexes rather than full scans.
//
// Run with:
//   docker exec <container> cypher-shell -u neo4j -p password < kg/schema.cypher
// ============================================================

// --- Scaffold node types (exported from Reactome) ---

CREATE CONSTRAINT rgp_identifier IF NOT EXISTS
    FOR (n:ReferenceGeneProduct) REQUIRE n.identifier IS UNIQUE;

CREATE CONSTRAINT pathway_stid IF NOT EXISTS
    FOR (n:Pathway) REQUIRE n.stId IS UNIQUE;

CREATE CONSTRAINT go_bp_identifier IF NOT EXISTS
    FOR (n:GO_BiologicalProcess) REQUIRE n.identifier IS UNIQUE;

CREATE CONSTRAINT ref_therapeutic_identifier IF NOT EXISTS
    FOR (n:ReferenceTherapeutic) REQUIRE n.identifier IS UNIQUE;

CREATE CONSTRAINT species_taxid IF NOT EXISTS
    FOR (n:Species) REQUIRE n.taxId IS UNIQUE;

CREATE CONSTRAINT compartment_identifier IF NOT EXISTS
    FOR (n:Compartment) REQUIRE n.identifier IS UNIQUE;

// --- Future node types (populated by the ChatEMT pipeline) ---

// InteractionEvent: one unique literature-extracted triplet (source, relationship_type, target)
//   Merged across all papers; evidence provenance lives in linked Evidence nodes.
//   Properties:
//     identifier             (str)  — triplet fingerprint (SHA256 of source|rel|target, 16 hex)
//     triplet_fingerprint    (str)  — same as identifier
//     relationship_type      (str)  — e.g. 'positive_causal', 'negative_causal', 'physical_interaction'
//     displayName            (str)  — human-readable summary
//     contradicting_evidence (str)  — 'true' if involved in a confirmed contradiction, else 'false'
CREATE CONSTRAINT interaction_event_identifier IF NOT EXISTS
    FOR (n:InteractionEvent) REQUIRE n.identifier IS UNIQUE;

// Evidence: one occurrence of a triplet extracted from a specific paper and text chunk
//   Properties:
//     identifier    (str) — stable ID (SHA256 of triplet_fingerprint|pmcid|chunk_doc_id, 16 hex)
//     evidence_text (str) — verbatim extracted text snippet
//     text_chunk_id (str) — MongoDB chunk_doc _id for the source text chunk
//     pmcid         (str) — PMC paper identifier
CREATE CONSTRAINT evidence_identifier IF NOT EXISTS
    FOR (n:Evidence) REQUIRE n.identifier IS UNIQUE;

// ExtractedTerm: novel term from NED that did not match Reactome or GO
//   Properties:
//     displayName (str) — term text
//     run_id      (str) — pipeline run identifier
//     cluster_id  (int) — refined_cluster integer from clustering pipeline
CREATE CONSTRAINT extracted_term_name IF NOT EXISTS
    FOR (n:ExtractedTerm) REQUIRE n.displayName IS UNIQUE;

CREATE CONSTRAINT extracted_term_cluster_id IF NOT EXISTS
    FOR (n:ExtractedTerm) REQUIRE n.cluster_id IS UNIQUE;

// Molecule: small molecule
//   Properties:
//     identifier  (str)
//     displayName (str)
CREATE CONSTRAINT molecule_identifier IF NOT EXISTS
    FOR (n:Molecule) REQUIRE n.identifier IS UNIQUE;

// --- Indexes on commonly queried non-unique properties ---

CREATE INDEX pathway_display_name IF NOT EXISTS
    FOR (n:Pathway) ON (n.displayName);

CREATE INDEX rgp_display_name IF NOT EXISTS
    FOR (n:ReferenceGeneProduct) ON (n.displayName);

CREATE INDEX go_bp_display_name IF NOT EXISTS
    FOR (n:GO_BiologicalProcess) ON (n.displayName);

CREATE INDEX compartment_display_name IF NOT EXISTS
    FOR (n:Compartment) ON (n.displayName);

CREATE INDEX interaction_event_rel_type IF NOT EXISTS
    FOR (n:InteractionEvent) ON (n.relationship_type);

CREATE INDEX evidence_pmcid IF NOT EXISTS
    FOR (n:Evidence) ON (n.pmcid);

// Note: geneName on ReferenceGeneProduct is a list property. Neo4j 4.4 list
// indexes support full-list equality only, not element membership. Queries
// using `'GENE' IN n.geneName` will perform a full scan on this property.

// --- Literature relationship types (populated by kg/build_literature_kg.py) ---
//
// subject       — (InteractionEvent)-[:subject]->(entity node)   [was: hasSource]
// object        — (InteractionEvent)-[:object]->(entity node)    [was: hasTarget]
// hasEvidence   — (InteractionEvent)-[:hasEvidence]->(Evidence)
// extractedGene — (ExtractedTerm)-[:extractedGene]->(ReferenceGeneProduct)
//
// Note: ExtractedTerm nodes are created ONLY for novel clusters (no Reactome/GO match).
// Clusters matched to a Reactome Pathway or GO_BiologicalProcess node resolve
// directly to those scaffold nodes — no ExtractedTerm intermediate.
