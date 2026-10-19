# EMT-InMap Companion Website

Interactive Streamlit companion application for **EMT-InMap**: time-resolved native complexome mapping and co-fractionation mass spectrometry (CF/MS) analysis of epithelial-mesenchymal transition (EMT) protein complexes.

## Features
- **Protein Complex & Assembly Explorer:** Browse regulated MCF10A ClusterOne protein assemblies across the EMT time course (Days 0–8).
- **PPI & Causal Networks:** Explore protein-protein interaction trajectories, binned signal-to-noise metrics, and integrated literature/pathway causal networks.
- **Structural Atlas & Phospho-Interfaces:** Interactive structural inspection and AlphaFold/PDB interface model mappings.
- **LLM Hypothesis Generator:** Optional AI-powered hypothesis generation and functional enrichment analysis over pre-computed tables and knowledge graph tools.

## Local Development & Running
1. Install dependencies:
   ```bash
   pip install -r 05_companion_website/requirements.txt
   ```
2. Run the Streamlit app:
   ```bash
   streamlit run 05_companion_website/streamlit_app.py
   ```
