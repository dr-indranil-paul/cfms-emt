import unittest
from unittest.mock import patch

import kg_tools


class ProviderTests(unittest.TestCase):
    def test_openai_is_available_in_provider_catalogue(self):
        self.assertIn("OpenAI", kg_tools.PROVIDER_MODELS)
        self.assertEqual(kg_tools.PROVIDER_ENV_KEY["OpenAI"], "OPENAI_API_KEY")
        models = kg_tools.PROVIDER_MODELS["OpenAI"]
        self.assertIn("gpt-5.4", models)
        self.assertIn("gpt-5.5", models)
        self.assertIn("gpt-5.6-terra", models)
        self.assertIn("gpt-5.6-sol", models)
        self.assertIn("gpt-5.6-luna", models)

    @patch("kg_tools.build_context", return_value="assembly context")
    @patch("kg_tools._call_openai_compatible", return_value="hypothesis")
    def test_openai_uses_direct_api_endpoint(self, call_api, _build_context):
        result, context = kg_tools.generate_hypothesis(
            "How does MON1B regulate EMT?",
            api_key="test-key",
            model="gpt-5.4",
            provider="OpenAI",
        )

        self.assertEqual(result, "hypothesis")
        self.assertEqual(context, "assembly context")
        self.assertEqual(
            call_api.call_args.kwargs["base_url"], "https://api.openai.com/v1"
        )
        self.assertEqual(call_api.call_args.args[:2], ("test-key", "gpt-5.4"))


if __name__ == "__main__":
    unittest.main()
