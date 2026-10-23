"""Baseline 2 with ChatGPT (OpenAI GPT-6 Astra via OpenRouter): same batches and prompt as the Claude sub-agent run.

The model gets the batch text as the user message; no tools, no web. Output JSON is saved per batch as
data/baseline2/out/gpt_<variant>_batchNN.json and scored by `src/21_baseline2_llm.py score`.
usage: OPENROUTER_KEY_FILE=... python src/21b_baseline2_openrouter.py [--model openai/gpt-6-astra]
"""
import json
import os
import re
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
B2DIR = ROOT / "data/baseline2"
URL = "https://openrouter.ai/api/v1/chat/completions"
MODEL = sys.argv[sys.argv.index("--model") + 1] if "--model" in sys.argv else "openai/gpt-6-astra"


def call(path):
    key = Path(os.environ["OPENROUTER_KEY_FILE"]).read_text().strip()
    out = B2DIR / "out" / f"gpt_{path.stem}.json"
    if out.exists():
        return path.name, "cached"
    body = {"model": MODEL, "temperature": 0, "response_format": {"type": "json_object"},
            "messages": [{"role": "system", "content": "You are an expert microbial physiologist. Answer only with "
                                                       "the requested JSON object. Do not use tools or outside lookups."},
                         {"role": "user", "content": path.read_text()}]}
    for i in range(4):
        try:
            req = urllib.request.Request(URL, data=json.dumps(body).encode(),
                                         headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=900) as r:
                resp = json.load(r)
            text = resp["choices"][0]["message"]["content"]
            text = re.sub(r"^```(json)?|```$", "", text.strip()).strip()
            data = json.loads(text)
            out.write_text(json.dumps(data))
            usage = resp.get("usage", {})
            (B2DIR / "out" / f"gpt_{path.stem}.usage.json").write_text(json.dumps({"model": resp.get("model"),
                                                                                    "usage": usage}))
            return path.name, f"ok {len(data)} genomes, tokens {usage.get('prompt_tokens')}/{usage.get('completion_tokens')}"
        except Exception as e:
            err = e
            time.sleep(20 * (i + 1))
    return path.name, f"error {err}"


def main():
    batches = sorted((B2DIR / "batches").glob("*_batch*.txt"))
    with ThreadPoolExecutor(4) as ex:
        for name, status in ex.map(call, batches):
            print(name, status, flush=True)


if __name__ == "__main__":
    main()
