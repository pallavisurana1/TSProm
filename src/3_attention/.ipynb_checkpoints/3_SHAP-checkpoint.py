import os, json, numpy as np, pandas as pd, shap, torch
from transformers import AutoTokenizer, AutoModelForSequenceClassification

# ──────────────────────────── Helpers ────────────────────────────
def set_seed(seed: int = 42):
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

def load_job_from_json(path, index=0):
    with open(path, "r") as f:
        obj = json.load(f)
    if isinstance(obj, list):
        job = obj[index]
    elif isinstance(obj, dict):
        job = obj
    else:
        raise ValueError("jobs.json must be dict or list")
    return job


def load_sequences(input_csv=None, seqs_txt=None, sequence_col="sequence", top_n=None):
    """
    If top_n is None → returns ALL sequences.
    Otherwise returns top-N by |prob_1 - prob_0|.
    """
    if seqs_txt and os.path.exists(seqs_txt):
        with open(seqs_txt) as f:
            return [ln.strip() for ln in f if ln.strip()]

    if input_csv and os.path.exists(input_csv):
        df = pd.read_csv(input_csv)
        if sequence_col not in df.columns:
            raise ValueError(f"Column {sequence_col} not in {input_csv}")
        # find prob columns
        p0 = p1 = None
        for c in df.columns:
            lc = c.lower().strip()
            if lc in ("prob_0","p0","prob0"): p0 = c
            if lc in ("prob_1","p1","prob1"): p1 = c
        if p0 is not None and p1 is not None and top_n is not None:
            df["_confidence"] = (df[p1] - df[p0]).abs()
            df = df.sort_values("_confidence", ascending=False).head(top_n)
        return df[sequence_col].astype(str).tolist()

    raise ValueError("Provide either input_csv or seqs_txt")


def build_predict_fn_idsafe(tokenizer, model, device, target_class=1):
    mask_id = getattr(tokenizer, "mask_token_id", None)
    pad_id  = getattr(tokenizer, "pad_token_id", None)
    unk_id  = getattr(tokenizer, "unk_token_id", 0)
    if mask_id is None:
        mask_id = pad_id if pad_id is not None else unk_id

    def predict_prob(token_arrays):
        arr = np.array(token_arrays, dtype=object)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)

        encoded_examples = []
        for row in arr:
            tok_ids = []
            for t in list(row):
                if t is None or t == "":
                    tok_ids.append(unk_id)
                elif isinstance(t, str):
                    tid = tokenizer.convert_tokens_to_ids(t)
                    tok_ids.append(unk_id if tid is None else tid)
                else:
                    tok_ids.append(int(t))  # already an id
            enc = tokenizer.prepare_for_model(tok_ids, add_special_tokens=True, truncation=True)
            encoded_examples.append(enc)

        batch = tokenizer.pad(encoded_examples, return_tensors="pt").to(device)
        with torch.no_grad():
            probs = torch.softmax(model(**batch).logits, dim=-1)[:, target_class]
        return probs.detach().cpu().numpy()

    mask_token = tokenizer.mask_token or tokenizer.pad_token or "[PAD]"
    return predict_prob, mask_token


def explain_sequence(seq, tokenizer, predict_fn, mask_token, max_len=512, nsamples=800):
    raw_tokens = tokenizer.tokenize(seq)[:max_len]
    if not raw_tokens:
        print("[WARN] Empty tokenization; skipping.")
        return None, None
    tokens = []
    for t in raw_tokens:
        if t is None or t == "":
            tokens.append(mask_token)
        elif isinstance(t, bytes):
            tokens.append(t.decode(errors="ignore"))
        else:
            tokens.append(str(t))

    background = np.array([[mask_token]*len(tokens)], dtype=object)
    explainer = shap.KernelExplainer(predict_fn, background, link="logit")
    shap_vals = explainer.shap_values(np.array(tokens, dtype=object), nsamples=nsamples)
    if isinstance(shap_vals, list):
        shap_vals = np.array(shap_vals[0])

    ev = explainer.expected_value
    if isinstance(ev, (list, np.ndarray)):
        ev = np.array(ev).reshape(-1)[0]

    ex = shap.Explanation(values=shap_vals, base_values=ev, data=tokens)
    df_tokens = pd.DataFrame({"token": tokens, "shap_value": shap_vals})
    return ex, df_tokens


jobs_path = "/data/private/psurana/TSpDNA2/src/2_predict/jobs.json"
job = load_job_from_json(jobs_path, index=7)
job


checkpoint = job["model_path"]
output_dir = os.path.join(job["res_pdir"], "shap")
os.makedirs(output_dir, exist_ok=True)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
tokenizer_name = "jaandoui/DNABERT2-AttentionExtracted"

tokenizer = AutoTokenizer.from_pretrained(checkpoint, use_fast=True, trust_remote_code=True)
tokenizer = AutoTokenizer.from_pretrained(tokenizer_name, trust_remote_code=True)
model = AutoModelForSequenceClassification.from_pretrained(checkpoint, trust_remote_code=True).to(device).eval()


df_res=pd.read_csv(job["res_pdir"] + '/4sig_motifs_prePWM.csv')
df_res.head(2)

seqs = load_sequences(
    input_csv=job["res_pdir"] + '/4sig_motifs_prePWM.csv',
    sequence_col="Sequence", top_n=None)
len(seqs)


predict_fn, mask_token = build_predict_fn_idsafe(tokenizer, model, device, target_class=1)

seq = seqs[1]
tok = tokenizer.tokenize(seq)[:10]; print(len(tok))
bg = np.array([[mask_token]*len(tok)], dtype=object)

ex, df_tokens = explain_sequence(seq, tokenizer, predict_fn, mask_token, max_len=10, nsamples=10)

shap.plots.waterfall(ex, max_display=5)  