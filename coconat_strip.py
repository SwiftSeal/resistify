import torch
import esm
from transformers import T5EncoderModel, T5Tokenizer
import regex as re

sequence_id = 'zar1'
sequence = 'MVDAVVTVFLEKTLNILEEKGRTVSDYRKQLEDLQSELKYMQSFLKDAERQKRTNETLRTLVADLRELVYEAEDILVDCQLADGDDGNEQRSSNAWLSRLHPARVPLQYKKSKRLQEINERITKIKSQVEPYFEFITPSNVGRDNGTDRWSSPVYDHTQVVGLEGDKRKIKEWLFRSNDSQLLIMAFVGMGGLGKTTIAQEVFNDKEIEHRFERRIWVSVSQTFTEEQIMRSILRNLGDASVGDDIGTLLRKIQQYLLGKRYLIVMDDVWDKNLSWWDKIYQGLPRGQGGSVIVTTRSESVAKRVQARDDKTHRPELLSPDNSWLLFCNVAFAANDGTCERPELEDVGKEIVTKCKGLPLTIKAVGGLLLCKDHVYHEWRRIAEHFQDELRGNTSETDNVMSSLQLSYDELPSHLKSCILTLSLYPEDCVIPKQQLVHGWIGEGFVMWRNGRSATESGEDCFSGLTNRCLIEVVDKTYSGTIITCKIHDMVRDLVIDIAKKDSFSNPEGLNCRHLGISGNFDEKQIKVNHKLRGVVSTTKTGEVNKLNSDLAKKFTDCKYLRVLDISKSIFDAPLSEILDEIASLQHLACLSLSNTHPLIQFPRSMEDLHNLQILDASYCQNLKQLQPCIVLFKKLLVLDMTNCGSLECFPKGIGSLVKLEVLLGFKPARSNNGCKLSEVKNLTNLRKLGLSLTRGDQIEEEELDSLINLSKLMSISINCYDSYGDDLITKIDALTPPHQLHELSLQFYPGKSSPSWLSPHKLPMLRYMSICSGNLVKMQEPFWGNENTHWRIEGLMLSSLSDLDMDWEVLQQSMPYLRTVTANWCPELESFAIEDVGFRGGVWMKTPLHRT'

print("Loading ProtT5 model")
model = T5EncoderModel.from_pretrained("coconat-plms/prot_t5_xl_uniref50", weights_only=False)
tokenizer = T5Tokenizer.from_pretrained("coconat-plms/prot_t5_xl_uniref50", weights_only=False)
print("Done")

print("Embedding sequence")
sequence = re.sub(r"[UZOB]", "X", sequence)
ids = tokenizer.batch_encode_plus(sequence, add_special_tokens=True, padding="longest")
input_ids = torch.tensor(ids['input_ids'])
attention_mask = torch.tensor(ids['attention_mask'])
with torch.no_grad():
    embedding_repr = model(input_ids=input_ids, attention_mask=attention_mask)

print("Extracting embeddings")
emb = embedding_repr.last_hidden_state[0, :len(sequence)]
ret = emb.detach().cpu().numpy()

#ESM2

print("Loading ESM2 model")
model, alphabet = esm.pretrained.load_model_and_alphabet("coconat-plms/esm2/esm2_t33_650M_UR50D.pt")
print("Done")

batch_converter = alphabet.get_batch_converter()
batch_labels, batch_strs, batch_tokens = batch_converter([(sequence_id, sequence)])
tokens_len = (batch_tokens != alphabet.padding_idx).sum(1).item()
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=False)

token_representations = results["representations"][33]
embedding = token_representations[0, 1 : tokens_len - 1].detach().cpu().numpy()