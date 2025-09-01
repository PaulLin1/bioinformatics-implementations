# Pipeline

This is the location of the pipeline for streaming in big data.

# Pipeline

Pipelines for streaming in data. Both pipelines are highly modular for different algorithm and hyperparameters.

Originally, I was going to have pipeline that could injest different files and file types. But the two datasets, Kaggle Sequences and UniProt Trembl, where so different that I decided each one should have there one. There differences are highlighted in the section below.

---

## kaggle_pipeline

This pipeline is for a Kaggle sequence dataset I found online. The dataset consists of sequences and their classes. I'd assume this is for any ML tasks you want to do with the dataset, which makes sense given I found it on Kaggle. Therefore, the algorithms were used as feature extractors that could be used for downstream ML algorithms like SVM. 

Compared to the UniProt datset, this one is incredibly small. I mainly used this dataset for prototyping the UniProt one. As you can see in the root folder's README, the CPU implementation of the algorithm actually performs better than the CUDA version. I tried to optimized the CUDA version for a smaller dataset, but didn't spend too much time on it because I figured the CPU implementations would be way faster on that scale.

---

## pipeline

This pipeline is for the UniProt Trembl dataset. The reason I named it plainly pipeline is because I hope to use it for other datsets. The UniProt Trembl dataset contains mostly FASTA files which seems to be a standard format for sequences dataset. If I were to download more, I would tweak this pipeline to accomdate for that.

As of now, this pipeline is just used as a stress test for each algorithm. The data doesn't really mean anything to me, as it is just a bunch of proteins. I take a randomly chosen sequence and just compare it to every one in the dataset.

---

## Future Work

There were two issues I had when working on the project. First, I wish there was a dataset that was a medium size. I did do some testing on portions of the UniProt Trembl dataset, but having a medium (like 20 GB) dataset, would've been nice. Also, I wish there was a more clear task. The algorithm's I chose are mainly for processing and used for downstream tasks. It was pretty uneventful whenever I finished a run because the data didn't really mean anything. I guess that's why I enjoy ML more.
