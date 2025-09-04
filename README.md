# InformationAccretion
Compute Information Accretion (Information Content) of ontology terms from a dataset.

Information accretion ($ia$), introduced in [[Clark and Radivojac, 2013]](https://academic.oup.com/bioinformatics/article-pdf/29/13/i53/18535367/btt228.pdf), is a measure of how much information is added to an ontology annotation by node $v$ given that its parents $\mathcal{Pa}(v)$ are already annotated. Specifically, 

$$ia(v) = \log_2 \frac{1}{\mathrm{Pr}(v|\mathcal{Pa}(v))} = \log_2 \frac{\mathrm{Pr}(\mathcal{Pa}(v))}{\mathrm{Pr}(\mathcal{Pa}(v)|v) \mathrm{Pr}(v)}$$

To calculate this probability, we use the observed annotations in a dataset as the empirical distribution on which to estimate the probabilities. Thus, $\mathrm{Pr}(v|\mathcal{Pa}(v))$ is computed as the number of examples (proteins) annotated with the term corresponding to node $v$ *and* the terms corresponding to nodes $\mathcal{Pa}(v)$ divided by the number of examples annotated with the terms corresponding to nodes $\mathcal{Pa}(v)$. Note that the annotation of the term corresponding to node $v$ necessarily implies the annotation of terms $\mathcal{Pa}(v)$ to guarantee a *consistent subgraph* (valid annotations). This means that  $\mathrm{Pr}(\mathcal{Pa}(v)|v)=1$ and the denominator can be computed simply as the count of examples annotated with the term corresponding to node $v$.

Counts are regularized by introducing an artificial record (protein) annotated with every terms and adding it to the counts. That is, we add 1 to the counts for all terms to avoid dividing by 0.

In this code, we take the Gene Ontology graph in .obo format and a set of annotations and compute $ia$ for every term in the ontology.

If you are not sure that the terms have already been propagater of the structure of the graph has changed since that propagation, you may wish to re-propagate the term and/or check that they are already propagated before computing $ia$. 

[Clark and Radivojac, 2013] Clark WT, Radivojac P. Information-theoretic evaluation of predicted ontological annotations. Bioinformatics (2013) 29(13): i53-i61.

## Installation
Use pip to install:

```bash
pip install .
```

## Usage
After installation, the command `ia` is available from the CLI.

If you have a dataset of GO annotations, eg terms.tsv, and the corresponding GO graph structure, eg go-basic.obo, compute IA with the following call
```bash
ia --annot terms.tsv --graph go-basic.obo
```
This operation takes about 2 minutes for a dataset of 5M terms

If your annotated terms are not already propagated, enable the --prop flag:  
```
ia --annot terms.tsv --graph go-basic.obo --prop
```
This operation may add about 4 minutes to the runtime for a dataset of 5M terms

If `--graph` is omitted, the Gene Ontology graph is downloaded automatically.


The annotation file should be a tab-delimited file with headers `EntryID`, `term`, and `aspect`. 
```
EntryID term    aspect
target1 GO:0033058      BPO
target1 GO:0043056      BPO
target1 GO:0050879      BPO
target1 GO:0071965      BPO
target2 GO:0031987      BPO
target2 GO:0032501      BPO
```

The aspect should be one of BPO, CCO, or MFO. If the aspect column is missing, the code will look up the aspect. 
In this case, the file should include `EntryID` and `term` columns:
```
EntryID term
target1 GO:0033058
target1 GO:0043056
target1 GO:0050879
target1 GO:0071965
target2 GO:0031987
target2 GO:0032501
```


## Disclaimer
This is an independent implementation of the methods described in [Clark and Radivojac, 2013]

Please cite the orignal work [Clark and Radivojac, 2013] if using this measure and this repo if using this code.
