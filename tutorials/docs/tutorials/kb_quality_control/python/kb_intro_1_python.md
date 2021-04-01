<a href="https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/tutorials/docs/tutorials/kb_quality_control/python/kb_intro_1_python.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# Introduction to single-cell RNA-seq I: pre-processing and quality control

This Python notebook demonstrates the use of the kallisto and bustools programs for pre-processing single-cell RNA-seq data ([also available as an R notebook](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_intro_1_R.ipynb)). It streams in 1 million *C. elegans* reads, pseudoaligns them, and produces a *cells x genes* count matrix in about a minute. The notebook then performs some basic QC. It expands on a notebook prepared by Sina Booeshaghi for the Genome Informatics 2019 meeting, where he ran it in under 60 seconds during a 1 minute "lightning talk".

The [kallistobus.tools tutorials](https://www.kallistobus.tools/tutorials) site has a extensive list of follow-up tutorials and vignettes on single-cell RNA-seq.


```python
#@title 
from IPython.display import HTML
HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/x-rNofr88BM" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')

```




<iframe width="560" height="315" src="https://www.youtube.com/embed/x-rNofr88BM" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>



The notebook was written by A. Sina Booeshaghi and Lior Pachter. If you use the methods in this notebook for your analysis please cite the following publication, on which it is based:

* Melsted, P., Booeshaghi, A.S. et al. Modular and efficient pre-processing of single-cell RNA-seq. bioRxiv (2019). doi:10.1101/673285

##Setup


```python
# This is  used to time the running of the notebook
import time
start_time = time.time()
```

### Install python packages


```python
# These packages are pre-installed on Google Colab, but are included here to simplify running this notebook locally
%%capture
!pip install matplotlib
!pip install scikit-learn
!pip install numpy
!pip install scipy
```


```python
# Install packages for analysis and plotting
from scipy.io import mmread
from sklearn.decomposition import TruncatedSVD
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from scipy.sparse import csr_matrix
matplotlib.rcParams.update({'font.size': 22})
%config InlineBackend.figure_format = 'retina'
```


```python
%%time
%%capture
# `kb` is a wrapper for the kallisto and bustools program, and the kb-python package contains the kallisto and bustools executables.
!pip install kb-python==0.24.1
```

    CPU times: user 71.8 ms, sys: 20.8 ms, total: 92.6 ms
    Wall time: 9.89 s


### Download required files


```python
%%time
# The quantification of single-cell RNA-seq with kallisto requires an index. 
# Indices are species specific and can be generated or downloaded directly with `kb`. 
# Here we download a pre-made index for C. elegans (the idx.idx file) along with an auxillary file (t2g.txt) 
# that describes the relationship between transcripts and genes.
!wget -O idx.idx https://caltech.box.com/shared/static/82yv415pkbdixhzi55qac1htiaph9ng4.idx
!wget -O t2g.txt https://caltech.box.com/shared/static/cflxji16171skf3syzm8scoxkcvbl97x.txt
```

    --2021-03-31 23:35:57--  https://caltech.box.com/shared/static/82yv415pkbdixhzi55qac1htiaph9ng4.idx
    Resolving caltech.box.com (caltech.box.com)... 107.152.27.197
    Connecting to caltech.box.com (caltech.box.com)|107.152.27.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/82yv415pkbdixhzi55qac1htiaph9ng4.idx [following]
    --2021-03-31 23:35:57--  https://caltech.box.com/public/static/82yv415pkbdixhzi55qac1htiaph9ng4.idx
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/82yv415pkbdixhzi55qac1htiaph9ng4.idx [following]
    --2021-03-31 23:35:57--  https://caltech.app.box.com/public/static/82yv415pkbdixhzi55qac1htiaph9ng4.idx
    Resolving caltech.app.box.com (caltech.app.box.com)... 107.152.27.201
    Connecting to caltech.app.box.com (caltech.app.box.com)|107.152.27.201|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!TOjZC9YauQsUJ1QtiQcyl9nECELYGh_7GOcHslcOkuFjCXsoUakhYeEEEaN9FG0stqv7Uk-S8CQ9zwZqkbZT8yKcMu4qp9rDiMk9mPGdyAAt7wlMt5fGe8-ayhgpbUaa7JOrIUmD42YWwC8RVxtbw0wov8BcRoj9EHbA5HH9j4XQGFOxgVSP91uupfMciKBleiGzWCa4lHF7tEjYNHCkMPCKETj25x3X9VvLI8Tjp5IDsu_qoNHqpv53-owbBaaShO-DZGZwcVkhwqmz4iygFE1cjWaQxLGa46o1z48Yl2OBGo9VixhvEnoRAzMb0No9D5DZay6TiAPixZOV17GPzE793QTK0DGeTjjDEdxgwrZs6f5EBIPrZymUG7rj1HPEnXiAvqEpeZFj9IUvU6WovhB-V0QbFD-E4eSHZbcWK9T7CPE1Lv3gC1bhuYNULThU-bD9eb-qghz8-5qAjwfuEV4TsO7wMtbIUpsGQiMwhaaiO95Pwe5__pJNay9ErF16gA3V3vTvwwk8jUQsL0E1LMO4hX1pYGZ1gC1MTlAUfN1qzaKIm7gdmVtTlMlTOw1h5Jk3hQOOITeOYQeehRwE4X4kyZcXMvzK_8N59TKs7NZP_9PxXmG1r_jxWgfEbXxdaAOUQaaPfvZ1l-AlUZcMmO7snBIRT2oah0fe6VvAnz3Pl3IQ6dX1VKAsd0XPcfANf3NCj0ukSaT7Lj3gJNd82ReIEnYXFjkmHmB-eY84uXqjupKoSUXtX4YddjePPNqm6HPgoJPRSZnBcOyfACOTONkq4ppy6nujoqJviKOM3FYBXo7ng6HXf-mNeQSiqBR2ooAzgNKy6B3emVDdCX9kzZu7P-74doan52Z0BTs6qHQmriBDqO53HSOTswDluiYzOOneHRq73wgBJFBHh9GPELx0x1OGv0Anh59qn9M8gNtebnbqnChR2UonBBUg0qmK-MFWrwk3wJIOrlSl0mEcRKWIQqM6xIbLSX6AvJ-0MnZ0jBtHaCjfpSbZluXYSTDqMg6BLCpcDEhWk8WlPszIF4pOXHTrWhBDUY9L7typSxvmbj5IKL_imHwGxlFX3w9fshr7Z9PP-T30X8h8pQTJwrtI9Qtl_M6chYS235sM1QjiKUk3BcaYjtwUrJ6KoR2UVNJqoDtq8MOp386KHhvyAsuRanj_O_4BwDPhSog9gBwvdF5DPe8XrfSw7zNX2OXERLLbE61wDrjigDSeboYbT-ZwF4dkKKXFw9do4iPu6avsuDhv8t157NqmtWOV-31hai_N8EjueCaROPzlOBUglX09yWewZgAj-kvEfhPQtdmjdKcTl-Qv6V1yyhiEm-XJy9euoi0Sja-20BjguinCmBvodDopu-DJRkPBxYIGchtBQURBkK1GjkBD97eXHCQlYY0p4zpn_Q../download [following]
    --2021-03-31 23:35:57--  https://public.boxcloud.com/d/1/b1!TOjZC9YauQsUJ1QtiQcyl9nECELYGh_7GOcHslcOkuFjCXsoUakhYeEEEaN9FG0stqv7Uk-S8CQ9zwZqkbZT8yKcMu4qp9rDiMk9mPGdyAAt7wlMt5fGe8-ayhgpbUaa7JOrIUmD42YWwC8RVxtbw0wov8BcRoj9EHbA5HH9j4XQGFOxgVSP91uupfMciKBleiGzWCa4lHF7tEjYNHCkMPCKETj25x3X9VvLI8Tjp5IDsu_qoNHqpv53-owbBaaShO-DZGZwcVkhwqmz4iygFE1cjWaQxLGa46o1z48Yl2OBGo9VixhvEnoRAzMb0No9D5DZay6TiAPixZOV17GPzE793QTK0DGeTjjDEdxgwrZs6f5EBIPrZymUG7rj1HPEnXiAvqEpeZFj9IUvU6WovhB-V0QbFD-E4eSHZbcWK9T7CPE1Lv3gC1bhuYNULThU-bD9eb-qghz8-5qAjwfuEV4TsO7wMtbIUpsGQiMwhaaiO95Pwe5__pJNay9ErF16gA3V3vTvwwk8jUQsL0E1LMO4hX1pYGZ1gC1MTlAUfN1qzaKIm7gdmVtTlMlTOw1h5Jk3hQOOITeOYQeehRwE4X4kyZcXMvzK_8N59TKs7NZP_9PxXmG1r_jxWgfEbXxdaAOUQaaPfvZ1l-AlUZcMmO7snBIRT2oah0fe6VvAnz3Pl3IQ6dX1VKAsd0XPcfANf3NCj0ukSaT7Lj3gJNd82ReIEnYXFjkmHmB-eY84uXqjupKoSUXtX4YddjePPNqm6HPgoJPRSZnBcOyfACOTONkq4ppy6nujoqJviKOM3FYBXo7ng6HXf-mNeQSiqBR2ooAzgNKy6B3emVDdCX9kzZu7P-74doan52Z0BTs6qHQmriBDqO53HSOTswDluiYzOOneHRq73wgBJFBHh9GPELx0x1OGv0Anh59qn9M8gNtebnbqnChR2UonBBUg0qmK-MFWrwk3wJIOrlSl0mEcRKWIQqM6xIbLSX6AvJ-0MnZ0jBtHaCjfpSbZluXYSTDqMg6BLCpcDEhWk8WlPszIF4pOXHTrWhBDUY9L7typSxvmbj5IKL_imHwGxlFX3w9fshr7Z9PP-T30X8h8pQTJwrtI9Qtl_M6chYS235sM1QjiKUk3BcaYjtwUrJ6KoR2UVNJqoDtq8MOp386KHhvyAsuRanj_O_4BwDPhSog9gBwvdF5DPe8XrfSw7zNX2OXERLLbE61wDrjigDSeboYbT-ZwF4dkKKXFw9do4iPu6avsuDhv8t157NqmtWOV-31hai_N8EjueCaROPzlOBUglX09yWewZgAj-kvEfhPQtdmjdKcTl-Qv6V1yyhiEm-XJy9euoi0Sja-20BjguinCmBvodDopu-DJRkPBxYIGchtBQURBkK1GjkBD97eXHCQlYY0p4zpn_Q../download
    Resolving public.boxcloud.com (public.boxcloud.com)... 107.152.27.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|107.152.27.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 625579580 (597M) [application/octet-stream]
    Saving to: ‘idx.idx’
    
    idx.idx             100%[===================>] 596.60M  16.6MB/s    in 29s     
    
    2021-03-31 23:36:27 (20.7 MB/s) - ‘idx.idx’ saved [625579580/625579580]
    
    --2021-03-31 23:36:27--  https://caltech.box.com/shared/static/cflxji16171skf3syzm8scoxkcvbl97x.txt
    Resolving caltech.box.com (caltech.box.com)... 107.152.27.197
    Connecting to caltech.box.com (caltech.box.com)|107.152.27.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/cflxji16171skf3syzm8scoxkcvbl97x.txt [following]
    --2021-03-31 23:36:27--  https://caltech.box.com/public/static/cflxji16171skf3syzm8scoxkcvbl97x.txt
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/cflxji16171skf3syzm8scoxkcvbl97x.txt [following]
    --2021-03-31 23:36:27--  https://caltech.app.box.com/public/static/cflxji16171skf3syzm8scoxkcvbl97x.txt
    Resolving caltech.app.box.com (caltech.app.box.com)... 107.152.27.201
    Connecting to caltech.app.box.com (caltech.app.box.com)|107.152.27.201|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!gWbsrvnYAgYXxnB8WTrshtz_ZRR1ZUSrew8lXt473jS8-q72mE2fJN3i3uKWzqAO9s76XhRxjtXyRvHVa_zgHrCqbSzCTu0e7bB_YpGrmN_rSyEsHsFXKhZoSScPrTVeJsgIV_9ghOrIBIp3MJd44R83oTG9Ve3NqXm9-dqMtZw0RHK60ZrBt0JJW0vkTlBnqwUapw2QEKgGLgiJBYC0ia4KEwdVp_L3kUilD7Q_Iln6PxoZWsNVTj6j1pTpx_LfQFU-YPwVDwTB4HIapelPSYMVxkDmlrcxzxZpB9PoQFh0KkhP14XXhewMUzm_7RvmH9ElybOvT2KeHhd0X84Ppx-2w_IQOePXCFwspJPzutXE60JX5WtCSRtSQTxAjcmIppjJ-wn7s5JkhWda17KWwyZweoCHEuH3DJVTsH_tYtydKlnil0tGxiwvUkxpuds3zWBx19Z9bM3cQ5zKno1bZj0NkAXyxIvVPf0NftkSr-wYCTxW3mO3L4lvTrV9gCatnSXGtFqhOe1QRdZYsy34rnGcoa5ztcDbdf38mB2ordUH0ZOCo8ML11DfWssAmkVuAhk-9AI6gO-nCOPVXM54oDvPYLyFOmw185kPnMK4WZWCPMteL-qmACHF49MR7pvUxxCuiNBoQq43nZUYiMw1gB-kEjzJhSM4u1E1Gq4iQqmmp95FWIaXBELLoo7JIHqvKHkVRvNy3sgvMn_pWOfa8oN1mXTPoq8BmH5mgj6th2scwnnJJVAU8V_eWzfELyB71_gWT7tk16cum9d6hh3Ex8W7BGo1X8jEtK4i6LvvwYtng1ufaxBECUS_GecmuuOMJodxUPaj_cADwABylJxPwLVlyHAfxcCNBo8hrgHSE0t5m3juwzJb1qNO-UqHXUw_omUlnvNeN6gmPR-aGjLFzXpkPf45_i-diX9v8f8IVrOV5m-MCEVZLyRdyx5Zlte0-AYvt1IXOvwPYxaaQvo5RaCfnXu7javy1IP0LMU_Z_VFKKbCOFZC7Ote6vUTLQFHdM-I3reaPW-ZQAcqeF5v_UkefYvaRgR6KaY2SUCjzczOh7E_hM-aJn2fSMZBpWS7ywpXIPm3Q5oHBUKO9Q5Mbge6zUmbBo7Os6hzahJ5Y1IazhBpHfD1wR9aZY-Cl7C6cPvFN2-ZJsOuYI7lKmG4fa7UpXdIDjN7SUQWG_9oUtkkfXOlpokRkFbmrXD_fiUCbOkWFwhqxaSJMKyHG3__sHhv_esrnCPMi_i5Dpfk9bll251hO074910bq77Hap0svlqQcItAalRZ5jVfxB72EGOG3qYfNKzxfvU91t339Bhc-TIDc2XxI8wkp994UjcsB61MCcIQnxyhXni0rvfw_8KTV288JdgxkymGJBE./download [following]
    --2021-03-31 23:36:28--  https://public.boxcloud.com/d/1/b1!gWbsrvnYAgYXxnB8WTrshtz_ZRR1ZUSrew8lXt473jS8-q72mE2fJN3i3uKWzqAO9s76XhRxjtXyRvHVa_zgHrCqbSzCTu0e7bB_YpGrmN_rSyEsHsFXKhZoSScPrTVeJsgIV_9ghOrIBIp3MJd44R83oTG9Ve3NqXm9-dqMtZw0RHK60ZrBt0JJW0vkTlBnqwUapw2QEKgGLgiJBYC0ia4KEwdVp_L3kUilD7Q_Iln6PxoZWsNVTj6j1pTpx_LfQFU-YPwVDwTB4HIapelPSYMVxkDmlrcxzxZpB9PoQFh0KkhP14XXhewMUzm_7RvmH9ElybOvT2KeHhd0X84Ppx-2w_IQOePXCFwspJPzutXE60JX5WtCSRtSQTxAjcmIppjJ-wn7s5JkhWda17KWwyZweoCHEuH3DJVTsH_tYtydKlnil0tGxiwvUkxpuds3zWBx19Z9bM3cQ5zKno1bZj0NkAXyxIvVPf0NftkSr-wYCTxW3mO3L4lvTrV9gCatnSXGtFqhOe1QRdZYsy34rnGcoa5ztcDbdf38mB2ordUH0ZOCo8ML11DfWssAmkVuAhk-9AI6gO-nCOPVXM54oDvPYLyFOmw185kPnMK4WZWCPMteL-qmACHF49MR7pvUxxCuiNBoQq43nZUYiMw1gB-kEjzJhSM4u1E1Gq4iQqmmp95FWIaXBELLoo7JIHqvKHkVRvNy3sgvMn_pWOfa8oN1mXTPoq8BmH5mgj6th2scwnnJJVAU8V_eWzfELyB71_gWT7tk16cum9d6hh3Ex8W7BGo1X8jEtK4i6LvvwYtng1ufaxBECUS_GecmuuOMJodxUPaj_cADwABylJxPwLVlyHAfxcCNBo8hrgHSE0t5m3juwzJb1qNO-UqHXUw_omUlnvNeN6gmPR-aGjLFzXpkPf45_i-diX9v8f8IVrOV5m-MCEVZLyRdyx5Zlte0-AYvt1IXOvwPYxaaQvo5RaCfnXu7javy1IP0LMU_Z_VFKKbCOFZC7Ote6vUTLQFHdM-I3reaPW-ZQAcqeF5v_UkefYvaRgR6KaY2SUCjzczOh7E_hM-aJn2fSMZBpWS7ywpXIPm3Q5oHBUKO9Q5Mbge6zUmbBo7Os6hzahJ5Y1IazhBpHfD1wR9aZY-Cl7C6cPvFN2-ZJsOuYI7lKmG4fa7UpXdIDjN7SUQWG_9oUtkkfXOlpokRkFbmrXD_fiUCbOkWFwhqxaSJMKyHG3__sHhv_esrnCPMi_i5Dpfk9bll251hO074910bq77Hap0svlqQcItAalRZ5jVfxB72EGOG3qYfNKzxfvU91t339Bhc-TIDc2XxI8wkp994UjcsB61MCcIQnxyhXni0rvfw_8KTV288JdgxkymGJBE./download
    Resolving public.boxcloud.com (public.boxcloud.com)... 107.152.27.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|107.152.27.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 1010392 (987K) [text/plain]
    Saving to: ‘t2g.txt’
    
    t2g.txt             100%[===================>] 986.71K  4.08MB/s    in 0.2s    
    
    2021-03-31 23:36:30 (4.08 MB/s) - ‘t2g.txt’ saved [1010392/1010392]
    
    CPU times: user 302 ms, sys: 65.9 ms, total: 368 ms
    Wall time: 33.1 s


## Pseudoalignment and counting

In this notebook we pseudoalign 1 million *C. elegans* reads and count UMIs to produce a *cells x genes* matrix. Instead of being downloaded, the reads are streamed directly to the Google Colab notebook for quantification. 
See [this blog post](https://sinabooeshaghi.com/2019/07/09/fasterq-to-count-matrices-for-single-cell-rna-seq/) for more details on how the streaming works.

The data consists of a subset of reads from [GSE126954](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126954) described in the paper:

* Packer, J., Zhu, Q. et al. [A lineage-resolved molecular atlas of C. elegans embryogenesis at single-cell resolution](https://science.sciencemag.org/content/365/6459/eaax1971/tab-e-letters). Science (2019). doi:10.1126/science.aax1971

### Run kallisto and bustools


```python
%%time
# This step runs `kb` to quantify the reads. `kb` can take as input URLs where the reads are located, and will stream the data 
# to Google Colab where it is quantified as it is downloaded. This allows for quantifying very large datasets without first 
# downloading them and saving them to disk. 

!kb count -i idx.idx -g t2g.txt --overwrite -t 2 -x 10xv2 https://caltech.box.com/shared/static/fh81mkceb8ydwma3tlrqfgq22z4kc4nt.gz https://caltech.box.com/shared/static/ycxkluj5my7g3wiwhyq3vhv71mw5gmj5.gz
```

    [2021-03-31 23:36:30,948]    INFO Piping https://caltech.box.com/shared/static/fh81mkceb8ydwma3tlrqfgq22z4kc4nt.gz to tmp/fh81mkceb8ydwma3tlrqfgq22z4kc4nt.gz
    [2021-03-31 23:36:30,951]    INFO Piping https://caltech.box.com/shared/static/ycxkluj5my7g3wiwhyq3vhv71mw5gmj5.gz to tmp/ycxkluj5my7g3wiwhyq3vhv71mw5gmj5.gz
    [2021-03-31 23:36:30,954]    INFO Generating BUS file from
    [2021-03-31 23:36:30,954]    INFO         tmp/fh81mkceb8ydwma3tlrqfgq22z4kc4nt.gz
    [2021-03-31 23:36:30,954]    INFO         tmp/ycxkluj5my7g3wiwhyq3vhv71mw5gmj5.gz
    [2021-03-31 23:36:46,788]    INFO Sorting BUS file ./output.bus to tmp/output.s.bus
    [2021-03-31 23:36:49,288]    INFO Whitelist not provided
    [2021-03-31 23:36:49,288]    INFO Copying pre-packaged 10XV2 whitelist to .
    [2021-03-31 23:36:49,382]    INFO Inspecting BUS file tmp/output.s.bus
    [2021-03-31 23:36:49,827]    INFO Correcting BUS records in tmp/output.s.bus to tmp/output.s.c.bus with whitelist ./10xv2_whitelist.txt
    [2021-03-31 23:37:04,223]    INFO Sorting BUS file tmp/output.s.c.bus to ./output.unfiltered.bus
    [2021-03-31 23:37:06,631]    INFO Generating count matrix ./counts_unfiltered/cells_x_genes from BUS file ./output.unfiltered.bus
    CPU times: user 194 ms, sys: 30 ms, total: 224 ms
    Wall time: 37 s


### Exercises

- `kb` can quantify data that is streamed from a URL as in the example above, or can read in data from disk. Is it faster to stream data, or to download it first and then quantify it from disk?


```python
# %%time
# !wget https://caltech.box.com/shared/static/fh81mkceb8ydwma3tlrqfgq22z4kc4nt.gz 
# !wget https://caltech.box.com/shared/static/ycxkluj5my7g3wiwhyq3vhv71mw5gmj5.gz
# !kb count -i idx.idx -g t2g.txt --overwrite -t 2 -x 10xv2 fh81mkceb8ydwma3tlrqfgq22z4kc4nt.gz ycxkluj5my7g3wiwhyq3vhv71mw5gmj5.gz
```

- The -t option in `kb` sets the numnber of threads to be used. The Google Colab machine you are running on has two threads. If you run this notebook locally you can increase the number of threads beyond 2. As the number of threads is increased the running time decreases proportionately, although eventually the speed at which reads can be loaded from disk is a limiting factor. Verify that running `kb` with 1 thread on Google Colab takes about twice as long as with 2 threads.


```python
# %%time
# !kb count -i idx.idx -g t2g.txt --overwrite -t 1 -x 10xv2 https://caltech.box.com/shared/static/fh81mkceb8ydwma3tlrqfgq22z4kc4nt.gz https://caltech.box.com/shared/static/ycxkluj5my7g3wiwhyq3vhv71mw5gmj5.gz
```

## Basic QC 

### Represent the cells in 2D


```python
# Read in the count matrix that was output by `kb`.
mtx = mmread("/content/counts_unfiltered/cells_x_genes.mtx")
```


```python
# Perform SVD
tsvd = TruncatedSVD(n_components=2)
tsvd.fit(mtx)
X = tsvd.transform(mtx)
```


```python
# Plot the cells in the 2D PCA projection
fig, ax = plt.subplots(figsize=(10, 7))

ax.scatter(X[:,0], X[:,1], alpha=0.5, c="green")

plt.axis('off')
plt.show()
```


![png](kb_intro_1_python_files/kb_intro_1_python_25_0.png)


While the PCA plot shows the overall structure of the data, a visualization highlighting the density of points reveals a large number of droplets represented in the lower left corner.



```python
# density display for PCA plot
from scipy.interpolate import interpn

def density_scatter( x , y, ax = None, sort = True, bins = 20, **kwargs )   :
    """
    Scatter plot colored by 2d histogram
    """
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins)
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    sc = ax.scatter( x, y, c=z, **kwargs )
    return sc

fig, ax = plt.subplots(figsize=(7,7))

x = X[:,0]
y = X[:,1]

sc = density_scatter(x, y, ax=ax, cmap="Greens")

fig.colorbar(sc, ax=ax)
plt.axis('off')

plt.show()
```


![png](kb_intro_1_python_files/kb_intro_1_python_27_0.png)


The following plot helps clarify the reason for the concentrated points in the lower-left corner of the PCA plot.


```python
# Create sparse matrix representation of the count matrix
mtx = csr_matrix(mtx)
```

### Test for library saturation


```python
# Create a plot showing genes detected as a function of UMI counts.
fig, ax = plt.subplots(figsize=(10, 7))

ax.scatter(np.asarray(mtx.sum(axis=1))[:,0], np.asarray(np.sum(mtx>0, axis=1))[:,0], color="green", alpha=0.01)
ax.set_xlabel("UMI Counts")
ax.set_ylabel("Genes Detected")
ax.set_xscale('log')
ax.set_yscale('log', nonposy='clip')

ax.set_xlim((0.5, 4500))
ax.set_ylim((0.5,2000))


plt.show()
```


![png](kb_intro_1_python_files/kb_intro_1_python_31_0.png)


Here we see that there are a large number of near empty droplets. A useful approach to filtering out such data is the "knee plot" shown below.

### Examine the knee plot

The "knee plot" was introduced in the Drop-seq paper: 
- Macosko et al., [Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets](https://www.cell.com/fulltext/S0092-8674(15)00549-8), 2015. DOI:10.1016/j.cell.2015.05.002

In this plot cells are ordered by the number of UMI counts associated to them (shown on the *x*-axis), and the fraction of droplets with at least that number of cells is shown on the *y*-axis:


```python
# Create the "knee plot"
knee = np.sort((np.array(mtx.sum(axis=1))).flatten())[::-1]
fig, ax = plt.subplots(figsize=(10, 7))

ax.loglog(knee, range(len(knee)),linewidth=5, color="g")

ax.set_xlabel("UMI Counts")
ax.set_ylabel("Set of Barcodes")

plt.grid(True, which="both")
plt.show()

```


![png](kb_intro_1_python_files/kb_intro_1_python_34_0.png)



```python
# An option is to filter the cells and genes by a threshold
# row_mask = np.asarray(mtx.sum(axis=1)>30).reshape(-1)
# col_mask = np.asarray(mtx.sum(axis=0)>0).reshape(-1)
# mtx_filtered = mtx[row_mask,:][:,col_mask]
```

### Exercises

- The "knee plot" is sometimes shown with the UMI counts on the y-axis instead of the x-axis, i.e. flipped and rotated 90 degrees. Make the flipped and rotated plot. Is there a reason to prefer one orientation over the other?


```python
# # Create the flipped and rotated "knee plot"
# knee = np.sort((np.array(mtx.sum(axis=1))).flatten())[::-1]
# fig, ax = plt.subplots(figsize=(10, 7))
# 
# ax.loglog(range(len(knee)), knee,linewidth=5, color="g")
# 
# ax.set_xlabel("Set of Barcodes")
# ax.set_ylabel("UMI Counts")
# 
# plt.grid(True, which="both")
# plt.show()

```

For more information on this exercise see [Rotating the knee (plot) and related yoga](https://liorpachter.wordpress.com/2019/06/24/rotating-the-knee-plot-and-related-yoga/).

- The PCA subspaces form a [*flag*](https://en.wikipedia.org/wiki/Flag_(linear_algebra). This means, for example, that regardless of the number of dimensions chosen for the PCA dimensionality reduction, the 2D subspace remains the same. Verify this empirically. 

- As you increase the number of dimensions for the PCA reduction, you can also view the relationship between different suspaces. Explore this by changing the subspace dimensions visualized.




```python
#@title Exploring PCA subspsaces { run: "auto", vertical-output: true, display-mode: "both" }
n_components =  2#@param {type:"integer"}
dimension_A =  1#@param {type:"integer"}
dimension_B =  2#@param {type:"integer"}
# Perform SVD
tsvd = TruncatedSVD(n_components)
tsvd.fit(mtx)
X = tsvd.transform(mtx)

# Plot the cells in the 2D PCA projection
fig, ax = plt.subplots(figsize=(10, 7))

ax.scatter(X[:,dimension_A-1], X[:,dimension_B-1], alpha=0.5, c="green")

plt.axis('off')
plt.show()
```


![png](kb_intro_1_python_files/kb_intro_1_python_41_0.png)


## Discussion

This notebook has demonstrated the pre-processing required for single-cell RNA-seq analysis. `kb` is used to pseudoalign reads and to generate a *cells x genes* matrix. Following generation of a matrix, basic QC helps to assess the quality of the data.


```python
# Running time of the notebook
print("{:.2f} seconds".format((time.time()-start_time)))
```

    96.95 seconds


**Feedback**: please report any issues, or submit pull requests for improvements, in the [Github repository where this notebook is located](https://github.com/pachterlab/kallistobustools/).


```python

```
