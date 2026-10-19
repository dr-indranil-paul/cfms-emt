"""
Various utility functions for manipulating the knowledge graph
"""

# import statements
import pandas as pd
import numpy as np
import ast


# function definitions
def expand_dict_column(df, column_name):
    def parse_cell(cell):
        if isinstance(cell, dict):
            return cell  # Already a dict
        elif isinstance(cell, str):
            try:
                return ast.literal_eval(cell)  # Try parsing string to dict
            except (ValueError, SyntaxError):
                return {}  # Return empty dict if parsing fails
        else:
            return {}  # If not str or dict, treat as empty dict

    parsed_dicts = df[column_name].apply(parse_cell)
    dict_df = pd.json_normalize(parsed_dicts)
    result = pd.concat([df.drop(columns=[column_name]), dict_df], axis=1)

    return result


def ensure_list(val):
    # Handle missing values
    if val is None:
        return []
    if isinstance(val, float) and pd.isna(val):
        return []

    # If it's already list-like
    if isinstance(val, (list, tuple, set)):
        return list(val)

    # If it's a numpy array
    if isinstance(val, np.ndarray):
        return val.tolist()

    # If it's a string, try to parse as a list
    if isinstance(val, str):
        try:
            parsed = ast.literal_eval(val)
            if isinstance(parsed, (list, tuple, set)):
                return list(parsed)
            else:
                return [parsed]
        except (ValueError, SyntaxError):
            return [val]  # fallback: single string

    # Fallback: wrap in a list
    return [val]


def add_assemblies(kg_nodes_df, kg_edges_df, ppi_nodes_df, ppi_edges_df):
    """


    Assembly class instances will have 'node_labels' column in nodes_df that looks
    like this:

    ['DatabaseObject', 'Event', 'Assembly']

    Interaction -- UndirectedInteraction




    For each assembly, what needs to be done?



    pre-compute existing EWAS and RGP

    if I need to add new genes to these, don't update that set; just add to genes_check and it will be understood that they're in now


    """

    kg_nodes_df["ewas_category"] = np.nan

    # DATA INFO
    # Current size
    # print("kg_nodes_df:", kg_nodes_df.shape)

    # --- ReferenceGeneProduct ---
    mask_rgp = kg_nodes_df["node_labels"].apply(
        lambda x: isinstance(x, (list, tuple, set)) and "ReferenceGeneProduct" in x
    )
    # print("RGP count:", mask_rgp.sum())              # number of matching rows
    rgp_df_all = kg_nodes_df[mask_rgp]

    # Filter RGP nodes: keep only those with exactly 4 node_labels for linking
    rgp_df = rgp_df_all[
        rgp_df_all["node_labels"].apply(lambda x: len(ensure_list(x)) == 4)
    ]

    # Count different types of RGP nodes
    rgp_nodes_4_labels = len(rgp_df)
    rgp_nodes_5_labels = len(
        rgp_df_all[rgp_df_all["node_labels"].apply(lambda x: len(ensure_list(x)) == 5)]
    )
    rgp_nodes_other = len(rgp_df_all) - rgp_nodes_4_labels - rgp_nodes_5_labels

    print(f"RGP nodes analysis:")
    print(f"  Total RGP nodes: {len(rgp_df_all)}")
    print(f"  RGP nodes with 4 labels (will be used for linking): {rgp_nodes_4_labels}")
    print(
        f"  RGP nodes with 5 labels (kept in KG, not used for linking): {rgp_nodes_5_labels}"
    )
    if rgp_nodes_other > 0:
        print(f"  RGP nodes with other label counts: {rgp_nodes_other}")

    # print("RGP df shape:", rgp_df.shape)             # shape of the filtered df

    # --- EntityWithAccessionedSequence ---
    mask_ewas = kg_nodes_df["node_labels"].apply(
        lambda x: (
            isinstance(x, (list, tuple, set)) and "EntityWithAccessionedSequence" in x
        )
    )
    # print("EWAS count:", mask_ewas.sum())
    ewas_df = kg_nodes_df[mask_ewas]
    # print("EWAS df shape:", ewas_df.shape)

    # print(rgp_df)
    # print(rgp_df['displayName'])
    # print(rgp_df['geneName'])

    # print(ewas_df)

    genes_check = set()
    rgps_check = set()

    unique_rgps = set(x for row in rgp_df["geneName"].map(ensure_list) for x in row)
    # print("Unique RGPs:", unique_rgps)

    # RGP Analysis - Check for duplicates and overlaps
    # print("\n=== RGP Analysis ===")
    # print(f"Total RGP nodes: {len(rgp_df)}")

    # Flatten all gene names from all RGP nodes
    all_rgp_genes = []
    rgp_gene_to_nodes = {}  # Track which nodes contain each gene

    for idx, row in rgp_df.iterrows():
        gene_list = ensure_list(row["geneName"])
        for gene in gene_list:
            all_rgp_genes.append(gene)
            if gene not in rgp_gene_to_nodes:
                rgp_gene_to_nodes[gene] = []
            rgp_gene_to_nodes[gene].append(row["node_id"])

    # print(f"Total gene names across all RGPs: {len(all_rgp_genes)}")
    # print(f"Unique gene names: {len(unique_rgps)}")

    # Find genes that appear in multiple RGP nodes
    duplicate_genes = {
        gene: nodes for gene, nodes in rgp_gene_to_nodes.items() if len(nodes) > 1
    }

    # print(f"Genes appearing in multiple RGP nodes: {len(duplicate_genes)}")
    if duplicate_genes:
        print("Examples of duplicate genes:")
        for i, (gene, node_ids) in enumerate(list(duplicate_genes.items())[:5]):
            print(
                f"  {gene}: appears in {len(node_ids)} RGP nodes (node_ids: {node_ids})"
            )
        if len(duplicate_genes) > 5:
            print(f"  ... and {len(duplicate_genes) - 5} more")

    # Check RGP nodes with multiple genes
    multi_gene_rgps = []
    for idx, row in rgp_df.iterrows():
        gene_list = ensure_list(row["geneName"])
        if len(gene_list) > 1:
            multi_gene_rgps.append((row["node_id"], gene_list))

    print(f"RGP nodes containing multiple genes: {len(multi_gene_rgps)}")
    if multi_gene_rgps:
        print("Examples of multi-gene RGPs:")
        for i, (node_id, genes) in enumerate(multi_gene_rgps[:5]):
            print(f"  Node {node_id}: {genes}")
        if len(multi_gene_rgps) > 5:
            print(f"  ... and {len(multi_gene_rgps) - 5} more")

    print("=== End RGP Analysis ===\n")

    # Analysis: Check overlap between duplicate RGP genes and PPI genes
    print("=== PPI-RGP Overlap Analysis ===")
    ppi_genes = set(ppi_nodes_df["gene"].unique())
    duplicate_gene_names = set(duplicate_genes.keys())

    # Find overlap between duplicate RGP genes and PPI genes
    duplicate_rgp_in_ppi = duplicate_gene_names.intersection(ppi_genes)

    print(f"Total genes in PPI data: {len(ppi_genes)}")
    print(f"Duplicate RGP genes: {len(duplicate_gene_names)}")
    print(f"Duplicate RGP genes that appear in PPI data: {len(duplicate_rgp_in_ppi)}")
    print(
        f"Percentage of duplicate RGP genes in PPI: {len(duplicate_rgp_in_ppi) / len(duplicate_gene_names) * 100:.1f}%"
    )

    if duplicate_rgp_in_ppi:
        print("Examples of duplicate RGP genes in PPI data:")
        for i, gene in enumerate(list(duplicate_rgp_in_ppi)[:10]):
            rgp_nodes = duplicate_genes[gene]
            print(f"  {gene}: {len(rgp_nodes)} RGP nodes (node_ids: {rgp_nodes})")
        if len(duplicate_rgp_in_ppi) > 10:
            print(f"  ... and {len(duplicate_rgp_in_ppi) - 10} more")

    print("=== End PPI-RGP Overlap Analysis ===\n")

    # check_gene = 'MICAL3'
    # print(check_gene in unique_rgps)

    # print(kg_nodes_df.columns)
    # print(kg_nodes_df.head(5))

    # unique node labels
    unique_node_labels = set(
        x for row in kg_nodes_df["node_labels"].map(ensure_list) for x in row
    )
    # print("Unique node labels:", unique_node_labels)
    # print('Interaction' in unique_node_labels)
    # print('InteractionEvent' in unique_node_labels)

    # check EWAS node labels:
    # a = kg_nodes_df[kg_nodes_df['node_labels'].apply(lambda x: isinstance(x, (list, tuple, set)) and 'EntityWithAccessionedSequence' in x)].head(1)['node_labels']
    # for x in a:
    #    print(x)

    """
    ['node_id', 'node_labels', 'name', 'neo4j', 'checksum', 'version',
       'dbId', 'isInDisease', 'speciesName', 'diagramWidth', 'hasDiagram',
       'diagramHeight', 'releaseDate', 'schemaClass', 'isInferred', 'hasEHLD',
       'stIdVersion', 'displayName', 'stId', 'doi', 'isChimeric',
       'endCoordinate', 'referenceType', 'startCoordinate', 'keyword',
       'isSequenceChanged', 'sequenceLength', 'secondaryIdentifier', 'url',
       'databaseName', 'chain', 'description', 'geneName', 'comment',
       'identifier', 'taxId', 'abbreviation', 'accessUrl', 'dateTime',
       'project', 'initial', 'surname', 'firstname', 'orcidId', 'address',
       'note', 'resourceIdentifier', 'definition', 'accession', 'synonym',
       'label', 'journal', 'title', 'pages', 'volume', 'year',
       'pubMedIdentifier', 'text', 'oldStId', 'trivial', 'formula',
       'minUnitCount', 'coordinate', 'ecNumber', 'releaseStatus',
       'variantIdentifier', 'maxUnitCount', 'secondCoordinate',
       'systematicName', 'isOrdered', 'inn', 'approved', 'approvalSource',
       'type', 'uniformResourceLocator', 'ISBN', 'chapterTitle',
       'otherIdentifier', 'endPositionInReferenceSequence',
       'startPositionInReferenceSequence', 'alteredAminoAcidFragment',
       'stoichiometryKnown']
    
    """

    """
    best order -- 

    create assembly -- done

    iterate through assembly_genes to ensure presence of correct EWAS and RGP nodes

    after iterating through assembly_genes to ensure presence of correct EWAS and RGP nodes, can select all 
    relevant edges from ppi_edges_df and create an InteractionEvent for each one. connect the correct EWAS 
    partners to the InteractionEvent node, and connect the Assembly to the InteractionEvent node. 

    Assembly -- hasEvent -- InteractionEvent -- partners -- EWAS -- referenceEntity -- RGP
    
    """

    # Get unique assemblies and genes
    assemblies = ppi_nodes_df["complex"].unique()

    print(len(assemblies), "unique assemblies found in PPI data")

    # Create assembly nodes in bulk
    assembly_rows = []
    start_node_id = kg_nodes_df.shape[0] + 1

    for i, assembly in enumerate(assemblies):
        assembly_row = pd.Series(index=kg_nodes_df.columns, dtype=object)
        assembly_row["node_labels"] = ["DatabaseObject", "Event", "Assembly"]
        assembly_row["name"] = assembly
        assembly_row["displayName"] = assembly
        assembly_row["node_id"] = start_node_id + i
        assembly_rows.append(assembly_row)

    # Add all assembly nodes at once
    if assembly_rows:
        assembly_df = pd.DataFrame(assembly_rows)
        kg_nodes_df = pd.concat([kg_nodes_df, assembly_df], ignore_index=True)

    # Get all unique genes that need EWAS nodes
    all_genes = ppi_nodes_df["gene"].unique()
    genes_needing_ewas = [gene for gene in all_genes if gene in unique_rgps]
    missing_rgps = [gene for gene in all_genes if gene not in unique_rgps]

    print(f"Total missing RGPs: {len(missing_rgps)}")
    if missing_rgps:
        print(
            f"Missing RGP genes: {missing_rgps[:10]}{'...' if len(missing_rgps) > 10 else ''}"
        )

    # Create EWAS nodes in bulk
    ewas_rows = []
    start_ewas_node_id = kg_nodes_df.shape[0] + 1

    for i, gene in enumerate(genes_needing_ewas):
        ewas_row = pd.Series(index=kg_nodes_df.columns, dtype=object)
        ewas_row["node_labels"] = [
            "DatabaseObject",
            "PhysicalEntity",
            "GenomeEncodedEntity",
            "EntityWithAccessionedSequence",
        ]
        ewas_row["name"] = gene
        ewas_row["displayName"] = gene
        ewas_row["node_id"] = start_ewas_node_id + i
        ewas_row["compartment"] = ""
        ewas_row["ewas_category"] = "ppi_ewas"
        ewas_rows.append(ewas_row)

    # Add all EWAS nodes at once
    if ewas_rows:
        ewas_df = pd.DataFrame(ewas_rows)
        kg_nodes_df = pd.concat([kg_nodes_df, ewas_df], ignore_index=True)

    # Update genes_check with all processed genes
    genes_check.update(genes_needing_ewas)

    # Create InteractionEvent nodes for unique gene-gene interactions
    # Get unique interactions (source-target pairs) across all assemblies
    unique_interactions = ppi_edges_df[["source", "target"]].drop_duplicates()

    # Filter to only include interactions where both genes have EWAS nodes
    valid_interactions = unique_interactions[
        (unique_interactions["source"].isin(genes_needing_ewas))
        & (unique_interactions["target"].isin(genes_needing_ewas))
    ]

    # Create InteractionEvent nodes in bulk
    interaction_rows = []
    start_interaction_node_id = kg_nodes_df.shape[0] + 1

    for i, (_, row) in enumerate(valid_interactions.iterrows()):
        source_gene = row["source"]
        target_gene = row["target"]

        interaction_row = pd.Series(index=kg_nodes_df.columns, dtype=object)
        interaction_row["node_labels"] = ["DatabaseObject", "Event", "InteractionEvent"]
        interaction_row["name"] = f"{source_gene}_{target_gene}_interaction"
        interaction_row["displayName"] = f"{source_gene} - {target_gene} interaction"
        interaction_row["node_id"] = start_interaction_node_id + i
        interaction_row["source_gene"] = source_gene
        interaction_row["target_gene"] = target_gene
        interaction_rows.append(interaction_row)

    # Add all InteractionEvent nodes at once
    if interaction_rows:
        interaction_df = pd.DataFrame(interaction_rows)
        kg_nodes_df = pd.concat([kg_nodes_df, interaction_df], ignore_index=True)

    print(f"Created {len(assemblies)} assembly nodes")
    print(f"Created {len(genes_needing_ewas)} EWAS nodes")
    print(f"Created {len(interaction_rows)} InteractionEvent nodes")

    # Link EWAS nodes to RGP nodes
    ppi_ewas_mask = kg_nodes_df["ewas_category"] == "ppi_ewas"
    ppi_ewas_nodes = kg_nodes_df[ppi_ewas_mask]

    # Create edges between EWAS and RGP nodes
    ewas_rgp_edges = []
    links_made = 0
    rgps_not_found = []

    for _, ewas_node in ppi_ewas_nodes.iterrows():
        gene_name = ewas_node["name"]

        # Find corresponding RGP node(s) for this gene
        rgp_matches = rgp_df[
            rgp_df["geneName"].apply(
                lambda x: gene_name in ensure_list(x) if pd.notna(x) else False
            )
        ]

        if len(rgp_matches) > 0:
            # Take the first match if multiple exist
            rgp_node = rgp_matches.iloc[0]
            if len(rgp_matches) > 1:
                # print(f"Warning: Multiple RGP matches found for gene {gene_name}. Using the first match with node_id {rgp_node['node_id']}.")
                pass

            # Create edge from EWAS to RGP (referenceEntity relationship)
            edge_row = pd.Series(index=kg_edges_df.columns, dtype=object)
            edge_row["source_id"] = ewas_node["node_id"]
            edge_row["target_id"] = rgp_node["node_id"]
            edge_row["relationship_type"] = "referenceEntity"
            edge_row["edge_id"] = len(kg_edges_df) + len(ewas_rgp_edges) + 1

            ewas_rgp_edges.append(edge_row)
            links_made += 1
        else:
            rgps_not_found.append(gene_name)

    # Add all EWAS-RGP edges at once
    if ewas_rgp_edges:
        ewas_rgp_df = pd.DataFrame(ewas_rgp_edges)
        kg_edges_df = pd.concat([kg_edges_df, ewas_rgp_df], ignore_index=True)

    print(f"Created {links_made} EWAS-RGP links")
    if rgps_not_found:
        print(
            f"RGPs not found for {len(rgps_not_found)} genes: {rgps_not_found[:10]}{'...' if len(rgps_not_found) > 10 else ''}"
        )

    # Connect InteractionEvent nodes to EWAS nodes (partner relationship)
    print("Creating InteractionEvent-EWAS partner relationships...")

    # Get InteractionEvent nodes
    interaction_mask = kg_nodes_df["node_labels"].apply(
        lambda x: isinstance(x, (list, tuple, set)) and "InteractionEvent" in x
    )
    interaction_nodes = kg_nodes_df[interaction_mask]

    # Build lookup dictionary: gene_name -> EWAS node_id (for ppi_ewas only)
    gene_to_ewas = {}
    for _, ewas_node in ppi_ewas_nodes.iterrows():
        gene_name = ewas_node["name"]
        gene_to_ewas[gene_name] = ewas_node["node_id"]

    # Create partner edges between InteractionEvent and EWAS nodes
    interaction_ewas_edges = []
    partner_links_made = 0

    for _, interaction_node in interaction_nodes.iterrows():
        source_gene = interaction_node["source_gene"]
        target_gene = interaction_node["target_gene"]
        interaction_node_id = interaction_node["node_id"]

        # Create edge from InteractionEvent to source EWAS (partner relationship)
        if source_gene in gene_to_ewas:
            edge_row = pd.Series(index=kg_edges_df.columns, dtype=object)
            edge_row["source_id"] = interaction_node_id
            edge_row["target_id"] = gene_to_ewas[source_gene]
            edge_row["relationship_type"] = "partner"
            edge_row["edge_id"] = (
                len(kg_edges_df) + len(ewas_rgp_edges) + len(interaction_ewas_edges) + 1
            )
            interaction_ewas_edges.append(edge_row)
            partner_links_made += 1

        # Create edge from InteractionEvent to target EWAS (partner relationship)
        if target_gene in gene_to_ewas:
            edge_row = pd.Series(index=kg_edges_df.columns, dtype=object)
            edge_row["source_id"] = interaction_node_id
            edge_row["target_id"] = gene_to_ewas[target_gene]
            edge_row["relationship_type"] = "partner"
            edge_row["edge_id"] = (
                len(kg_edges_df) + len(ewas_rgp_edges) + len(interaction_ewas_edges) + 1
            )
            interaction_ewas_edges.append(edge_row)
            partner_links_made += 1

    # Add all InteractionEvent-EWAS edges at once
    if interaction_ewas_edges:
        interaction_ewas_df = pd.DataFrame(interaction_ewas_edges)
        kg_edges_df = pd.concat([kg_edges_df, interaction_ewas_df], ignore_index=True)

    print(f"Created {partner_links_made} InteractionEvent-EWAS partner links")

    # Create edges between Assembly nodes and InteractionEvent nodes (hasEvent relationship)
    print("Creating Assembly-InteractionEvent hasEvent relationships...")

    # Get Assembly nodes
    assembly_mask = kg_nodes_df["node_labels"].apply(
        lambda x: isinstance(x, (list, tuple, set)) and "Assembly" in x
    )
    assembly_nodes = kg_nodes_df[assembly_mask]

    # Create hasEvent edges between Assembly and InteractionEvent nodes
    assembly_interaction_edges = []
    hasEvent_links_made = 0

    for _, assembly_node in assembly_nodes.iterrows():
        assembly_name = assembly_node["name"]
        assembly_node_id = assembly_node["node_id"]

        # Get all edges in this assembly from ppi_edges_df
        assembly_edges = ppi_edges_df[ppi_edges_df["complex"] == assembly_name]

        # For each edge in this assembly, find the corresponding InteractionEvent
        for _, ppi_edge in assembly_edges.iterrows():
            source_gene = ppi_edge["source"]
            target_gene = ppi_edge["target"]

            # Find the InteractionEvent node for this gene pair
            # InteractionEvent names are formatted as "source_target_interaction"
            interaction_name = f"{source_gene}_{target_gene}_interaction"

            # Find the InteractionEvent node with this name
            matching_interaction = interaction_nodes[
                interaction_nodes["name"] == interaction_name
            ]

            if len(matching_interaction) > 0:
                interaction_node_id = matching_interaction.iloc[0]["node_id"]

                # Create edge from Assembly to InteractionEvent (hasEvent relationship)
                edge_row = pd.Series(index=kg_edges_df.columns, dtype=object)
                edge_row["source_id"] = assembly_node_id
                edge_row["target_id"] = interaction_node_id
                edge_row["relationship_type"] = "hasEvent"
                edge_row["edge_id"] = (
                    len(kg_edges_df) + len(assembly_interaction_edges) + 1
                )

                assembly_interaction_edges.append(edge_row)
                hasEvent_links_made += 1

    # Add all Assembly-InteractionEvent edges at once
    if assembly_interaction_edges:
        assembly_interaction_df = pd.DataFrame(assembly_interaction_edges)
        kg_edges_df = pd.concat(
            [kg_edges_df, assembly_interaction_df], ignore_index=True
        )

    print(f"Created {hasEvent_links_made} Assembly-InteractionEvent hasEvent links")

    # print(kg_edges_df.columns)

    kg_edges_df["source_id"] = kg_edges_df["source_id"].astype(int)
    kg_edges_df["target_id"] = kg_edges_df["target_id"].astype(int)

    return kg_nodes_df, kg_edges_df

    """
    ppi_edges_df contains edges between genes. for each assembly in 'assemblies', 
    consider the subset of rows in ppi_edges_df where 'complex' matches the assembly. 
    each row contains a 'source' gene and a 'target' gene. for each row, we need to create 
    an 'InteractionEvent' node in kg_nodes_df and link it to the assembly node. node labels 
    should be: 

    Each InteractionEvent node should be linked to the correct EWAS nodes for the 
    source and target genes. there may be multiple EWAS nodes for a given gene, so use the EWAS 
    node where 

    ['ewas_category'] == 'ppi_ewas'




    I now need to create edges between Assembly nodes and InteractionEvent nodes (hasEvent relationship)




    
    """

    """
    FUNCTION PSEUDOCODE: 

    for each assembly: 
    
        add assembly class 

        iterate through genes in assembly 

            check if needed RGP and EWAS nodes exist for each gene 


    
    ADD: 

    add assembly node for each assembly 

    need to make sure there is exactly one RGP match; not more than one (this should be okay, but verify) 

    for each gene in an assembly: 

        if it has rgp and ewas, continue (this will only occur if protein has already come up in another assembly) 

        if it has rgp but not ewas yet (all of them at first, because I am adding new EWAS), create new EWAS

            also link this new EWAS to RGP 

        if it needs RGP and EWAS (530) -- decide -- search, or add RGP


        once all needed EWAS nodes exist for an assembly, pull list of interactions -- edges involved in assembly 

        create new interactions (assembly edges) for all -- interactions will link to EWAS nodes and also to the assembly 



        write pseudocode in docstring, then implement 


    then come up with systematic way to check KG integrity, PPI data add checks 

    then load back into graph db

    can get version complete, then iterate to polish (e.g. missing RGP nodes) 


    """

    return kg_nodes_df, kg_edges_df
