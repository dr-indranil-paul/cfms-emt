import unittest

from flow_diagrams import mermaid_to_plotly, parse_mermaid_flowchart


USER_EXAMPLE = """flowchart TD
    A[EMT Induction] --> B[C0197 INTAC MON1B]
    B -| C[C0457 Hypoxia Copper]
    C --> D[Cell Stress DAMPs]
    D --> E[T cell Priming]
    B --> F[RAB7 Activation]
    F --> G[Checkpoint Cytokine]
    G --> H[T cell Suppression]
    H -| E
"""


class FlowDiagramTests(unittest.TestCase):
    def test_repairs_biological_inhibition_syntax(self):
        nodes, edges = parse_mermaid_flowchart(USER_EXAMPLE)

        self.assertEqual(len(nodes), 8)
        self.assertEqual(len(edges), 8)
        self.assertEqual((edges[-1].source, edges[-1].target), ("H", "E"))
        self.assertEqual(sum(edge.kind == "inhibition" for edge in edges), 2)

    def test_parses_strict_prompt_syntax_and_edge_labels(self):
        text = """flowchart TD
        A[EMT induction] -->|activates| B[C0197 INTAC MON1B]
        B -.->|inhibits| C[T cell suppression]
        """
        nodes, edges = parse_mermaid_flowchart(text)

        self.assertEqual(nodes["B"], "C0197 INTAC MON1B")
        self.assertEqual(edges[0].label, "activates")
        self.assertEqual(edges[0].kind, "activation")
        self.assertEqual(edges[1].label, "inhibits")
        self.assertEqual(edges[1].kind, "inhibition")

    def test_joins_wrapped_node_labels(self):
        text = """flowchart TD
        A[EMT
          induction] --> B[Immune response]
        """
        nodes, edges = parse_mermaid_flowchart(text)

        self.assertEqual(nodes["A"], "EMT induction")
        self.assertEqual(len(edges), 1)

    def test_builds_styled_plotly_figure(self):
        figure = mermaid_to_plotly(USER_EXAMPLE)

        self.assertIsNotNone(figure)
        self.assertEqual(figure.layout.paper_bgcolor, "#ffffff")
        self.assertGreaterEqual(figure.layout.height, 430)
        node_count = len(parse_mermaid_flowchart(USER_EXAMPLE)[0])
        self.assertGreaterEqual(len(figure.layout.shapes), node_count * 3)


if __name__ == "__main__":
    unittest.main()
