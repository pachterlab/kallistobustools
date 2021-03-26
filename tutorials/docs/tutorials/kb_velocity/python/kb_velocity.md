# Pre-processing and RNA velocity analysis of single-cell RNA-seq data with kallisto|bustools.

In this notebook, we will perform pre-processing and RNA velocity analysis of human week 10 fetal forebrain dataset ([SRR6470906](https://www.ncbi.nlm.nih.gov/sra/?term=SRR6470906) and [SRR6470907](https://www.ncbi.nlm.nih.gov/sra/?term=SRR6470907)) from [La Manno et al., 2018](https://doi.org/10.1038/s41586-018-0414-6) using the **kallisto | bustools** workflow, implemented with a wrapper called `kb`. It was developed by Kyung Hoi (Joseph) Min and A. Sina Booeshaghi.

__Note:__ The human RNA velocity index used in this tutorial requires at least 16GB of RAM, which means you need to be on a GPU runtime on Google Colab. To change your runtime to a GPU instance, go to `Runtime` > `Change runtime type` > select `GPU` from the `Hardware accelerator` dropdown.

Downloading and processing the data takes a while (> 2 hours). If you would like to use checkpoint files, which reduce this to ~30 minutes, follow these steps:
1. Skip the *Download the data* step.
2. Follow through the pre-processing steps. Stop before *Generate RNA velocity count matrices*.
3. Uncomment the commented commands in the cells that run `kb count`.
4. Continue through the tutorial.


```
!date
```

    Thu Dec 12 20:34:50 UTC 2019


## Pre-processing

### Download the data

__Note:__ We use the `-O` option for `wget` to rename the files to easily identify them.


```
%%time
!wget https://caltech.box.com/shared/static/8w79k2ydhqigkb4cfhosc6k32zycizdq.txt -O checksums.txt
# SRR6470907
!wget https://caltech.box.com/shared/static/nvzqphhklk1yx938l6omursw7sr68y43.gz -O SRR6470906_S1_L001_R1_001.fastq.gz
!wget https://caltech.box.com/shared/static/63fh2xa5t82x7s74rqa0e2u2ur59y5ox.gz -O SRR6470906_S1_L001_R2_001.fastq.gz
!wget https://caltech.box.com/shared/static/zqi3durukillaw1pbns1kd1lowyfg5qk.gz -O SRR6470906_S1_L002_R1_001.fastq.gz
!wget https://caltech.box.com/shared/static/i56qojfz41ns1kw9z86sla0vawsch96t.gz -O SRR6470906_S1_L002_R2_001.fastq.gz
# SRR6470907
!wget https://caltech.box.com/shared/static/vrditbbk38tw3f61fwpg504vcc5x09ci.gz -O SRR6470907_S1_L001_R1_001.fastq.gz
!wget https://caltech.box.com/shared/static/8ud3otwztjeqlmjctbu1fw7hg3k56ejr.gz -O SRR6470907_S1_L001_R2_001.fastq.gz
!wget https://caltech.box.com/shared/static/ln14jjd4tz3hvgxf8zj2kmokof7f1nrf.gz -O SRR6470907_S1_L002_R1_001.fastq.gz
!wget https://caltech.box.com/shared/static/o5bwf9u2g7egi02by3e3hbvov8fgwbb3.gz -O SRR6470907_S1_L002_R2_001.fastq.gz
```

    --2019-12-12 16:45:17--  https://caltech.box.com/shared/static/8w79k2ydhqigkb4cfhosc6k32zycizdq.txt
    Resolving caltech.box.com (caltech.box.com)... 103.116.4.197
    Connecting to caltech.box.com (caltech.box.com)|103.116.4.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/8w79k2ydhqigkb4cfhosc6k32zycizdq.txt [following]
    --2019-12-12 16:45:17--  https://caltech.box.com/public/static/8w79k2ydhqigkb4cfhosc6k32zycizdq.txt
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/8w79k2ydhqigkb4cfhosc6k32zycizdq.txt [following]
    --2019-12-12 16:45:18--  https://caltech.app.box.com/public/static/8w79k2ydhqigkb4cfhosc6k32zycizdq.txt
    Resolving caltech.app.box.com (caltech.app.box.com)... 103.116.4.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|103.116.4.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!HROZcLQR4SR8wMrrjExg1KUpkObpwaMUNr5KN5DXu3Zn3FrctuyWJpB8m11X3wFIGVNxY8WZqjNGudP-Nw-ajzZIz-1kTO7dtxXAzVUs_8gVstlxRn0dh7tofcn9SEpe0L9Bjohw_8kCgdAEJCWUmqOqODDDO0MJ48xbHi8Xs8a_R6WfralHUXnJvZaj1JAOEuvr1U8Vtsl2y2AzzKIuUspRnN0vK1_lro9z34Ydrjhp7_6O71TKRtbNnREVwJUrj6P3PQ7IQBUmOr9t49TI9I8wLHpTGajKNDc5bweuFCjggl_gfDfTqGjeVk8YtGM3j8sl8qzblSahy8lWAhOLwCiGgwcsw1W-cPKMZssb2Ij5Bq5WCagQM_K5xRhWLPYd_-rn90p6kVHAX323pRRFcYyotJPAZ8UXNmmQYa6_CMOhcH5agYtE4jEb6vKNQUUezYHEGkv8ghoM9LFjd6Wzc30c66blugFTD0hpVQcPQdLiA4gquc_BrpKLsfQ8B6MPfndlhEIE2uZpXHtE2HIEsiBpLWbUTOXV2vK6xtUVsVFPeq9xQyGGe3_5TFarqqNzGIwO68uGfZ8RneCScuQVbnDym7zzZTCja7jcrHVP5Lj3gOcG2SeK8l4srpLZsfyBBIPrZkHaZ4Zmlan9PK7YEFVssJU5Ul-_k8whB67050b1-dS6aLZ2v9_rnH1h2sk_ymTKI_8a-PKAvECCvNy8AxABvk6UD4A_FEWjLYhonYYogAc73_2EvNT8gq9dSKIr3n5q2kj2vizM0frtoEHHSa1e89pLYsyGWOp5oeX9gGXMW9o0iyiqXVdSgEOrYL-PX4rETSJs4UCpupzHGanKhkIvBwmzyEPblYc3jykdp5z8wDP72Z0cCRTh7KPeo8yRe7D2Qpz_KUgPxh2uHqUdQ2F7HLfN1PIVoO9BK7bPYx7BKEmqahEy2ji7jGl1CBjFrg6SBQQmqr5mxiouiiCSLTecUT-w2s9XxDj1arBeGA0xX0T6gZn9zueU8HaYa8IV_l6PFjNZVe7Cd6kr9DRG4zRJ_n4K2tqNQouUvTCAvEc461DX39godP_dEP1XjS535ZYi8UmnFMAhsqDFXEJ9Ib9QEs1aoew7pDPhPRqwgPKnMm5TkNrzHXMt2Ng0HsoLcixfLWXG_oxzfDrZ-C9WlmN47Z0IoNrt9Wmd0iYDo_AU4KTfV_LXG7vdULCXMMpE9hulQBNNrmJCR5bC00RXm62Sw-D_oe1PcTfpWehBmlHMj6AWfo36pfHyeI43Dcxlu3Ryvc8JLCM5DpDla0DzjZMgVnr_sBUMvT7_QYXh1SNUyt-PDeKEqiXYV8pc6meDDWXSGTHEdxpkFvOqYV82A3C6BYcU1RO2aJY4rpxQO-y6/download [following]
    --2019-12-12 16:45:19--  https://public.boxcloud.com/d/1/b1!HROZcLQR4SR8wMrrjExg1KUpkObpwaMUNr5KN5DXu3Zn3FrctuyWJpB8m11X3wFIGVNxY8WZqjNGudP-Nw-ajzZIz-1kTO7dtxXAzVUs_8gVstlxRn0dh7tofcn9SEpe0L9Bjohw_8kCgdAEJCWUmqOqODDDO0MJ48xbHi8Xs8a_R6WfralHUXnJvZaj1JAOEuvr1U8Vtsl2y2AzzKIuUspRnN0vK1_lro9z34Ydrjhp7_6O71TKRtbNnREVwJUrj6P3PQ7IQBUmOr9t49TI9I8wLHpTGajKNDc5bweuFCjggl_gfDfTqGjeVk8YtGM3j8sl8qzblSahy8lWAhOLwCiGgwcsw1W-cPKMZssb2Ij5Bq5WCagQM_K5xRhWLPYd_-rn90p6kVHAX323pRRFcYyotJPAZ8UXNmmQYa6_CMOhcH5agYtE4jEb6vKNQUUezYHEGkv8ghoM9LFjd6Wzc30c66blugFTD0hpVQcPQdLiA4gquc_BrpKLsfQ8B6MPfndlhEIE2uZpXHtE2HIEsiBpLWbUTOXV2vK6xtUVsVFPeq9xQyGGe3_5TFarqqNzGIwO68uGfZ8RneCScuQVbnDym7zzZTCja7jcrHVP5Lj3gOcG2SeK8l4srpLZsfyBBIPrZkHaZ4Zmlan9PK7YEFVssJU5Ul-_k8whB67050b1-dS6aLZ2v9_rnH1h2sk_ymTKI_8a-PKAvECCvNy8AxABvk6UD4A_FEWjLYhonYYogAc73_2EvNT8gq9dSKIr3n5q2kj2vizM0frtoEHHSa1e89pLYsyGWOp5oeX9gGXMW9o0iyiqXVdSgEOrYL-PX4rETSJs4UCpupzHGanKhkIvBwmzyEPblYc3jykdp5z8wDP72Z0cCRTh7KPeo8yRe7D2Qpz_KUgPxh2uHqUdQ2F7HLfN1PIVoO9BK7bPYx7BKEmqahEy2ji7jGl1CBjFrg6SBQQmqr5mxiouiiCSLTecUT-w2s9XxDj1arBeGA0xX0T6gZn9zueU8HaYa8IV_l6PFjNZVe7Cd6kr9DRG4zRJ_n4K2tqNQouUvTCAvEc461DX39godP_dEP1XjS535ZYi8UmnFMAhsqDFXEJ9Ib9QEs1aoew7pDPhPRqwgPKnMm5TkNrzHXMt2Ng0HsoLcixfLWXG_oxzfDrZ-C9WlmN47Z0IoNrt9Wmd0iYDo_AU4KTfV_LXG7vdULCXMMpE9hulQBNNrmJCR5bC00RXm62Sw-D_oe1PcTfpWehBmlHMj6AWfo36pfHyeI43Dcxlu3Ryvc8JLCM5DpDla0DzjZMgVnr_sBUMvT7_QYXh1SNUyt-PDeKEqiXYV8pc6meDDWXSGTHEdxpkFvOqYV82A3C6BYcU1RO2aJY4rpxQO-y6/download
    Resolving public.boxcloud.com (public.boxcloud.com)... 103.116.4.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|103.116.4.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 552 [text/plain]
    Saving to: ‘checksums.txt’
    
    checksums.txt       100%[===================>]     552  --.-KB/s    in 0s      
    
    2019-12-12 16:45:19 (146 MB/s) - ‘checksums.txt’ saved [552/552]
    
    --2019-12-12 16:45:21--  https://caltech.box.com/shared/static/nvzqphhklk1yx938l6omursw7sr68y43.gz
    Resolving caltech.box.com (caltech.box.com)... 103.116.4.197
    Connecting to caltech.box.com (caltech.box.com)|103.116.4.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/nvzqphhklk1yx938l6omursw7sr68y43.gz [following]
    --2019-12-12 16:45:21--  https://caltech.box.com/public/static/nvzqphhklk1yx938l6omursw7sr68y43.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/nvzqphhklk1yx938l6omursw7sr68y43.gz [following]
    --2019-12-12 16:45:21--  https://caltech.app.box.com/public/static/nvzqphhklk1yx938l6omursw7sr68y43.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 103.116.4.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|103.116.4.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!jD0baoEwhAs7CayO95UyhSNszeFfA_EqNIM7jM-ZixrgAG6UaRFMLlmYPaJ8YVAXdFfxrrHnQc0CJHNUY_6RT3QrlkhFFuq92WCeuhmH1Ltvb75ni-GO-K2Ltbni0HMOij3IAWXLEz4IqZLNeMjY1Ln2JyW8cLN_63RoCKDAkqv8okgYR_sB92M6_v2_Cuj3UvfzH4u8H6DKzeFch8_SfF7gubKOAFFmgzG1X4z6tXQBAgf6ezIfkS-K20C3Uw5PK-wxFvUVLIHi_xm9W4FaubPBkxbVlBUqB0WFhMLOlKty3KcLKptS67fy8xMpYQGVmhAUCUzkmOK33rKOBSdePbM9XeicLsreHA1xq41eYCmE_S_3-7YOMFMoS-3wjEjsatR-yL_wlHR7Xop2I4Qr9-1JEpBye_C8TtUjLeC7OPzvxvnZJ4r4nWtY1qH9dTJODRf2B9twvXLoE8HE-MR3noULJvEArKEn5d2dIDk1wSpNnHPkKhOhNTp1RP8fjTQXBtUs8vbymUAhDX2lhq0pR4mjbR9Rwphvvom0hpeeuBFEI1R-Wn4M7GnAzdf0cqNwTf2OwuMIq5spmd2DAVRe7QHSSiHEkpd4BDlXg2JTuKIhX4l4c1naAM_3g1zD9EdBWN3pdHamjn3G5Gtk3KuIoMKuOy_fnvIo-nKnrO3DE4J0bgim-GdTQHRbsxBGu2MRU9AVVfgw_yaf2c6Hd_NkAd0GFFJgyEelWVfLYu214F0njqpRY6n3hqDRC_XRSRPbEH7-Xv0hEBYZzZqH31DKc3jb2sTPsG-jaF4Sv7WR8qOtjD2Xk6uv-1MWXefO-pw9GsLmgNGhWX85gmEbAHIxqKA_oJq9pSi3bFNbpWrimvxFGIEcr42b6Utg1oUP35pL5_T8_04SkIhxDihD9gxe6EoYmhpM2kzbcUyViEnG_RUBIFCj7kDaQb71ODQFYPxTXMyb3c_PipXe_GGm_JtZCp_xJeOOZZAs4QGbRXNaznpTzzyA3r66sSycCpbl18rkQQybOoHo8XNhkdzCt-cB3Dw_dkN_YJQ8_BkUIiiij3yETRPEO-NTvHNX_pqWWMRgxNb__RGzoMiANsGyqFlFZmQZnS03gEx6zwa7LhPsTxI6797c1ichD8uXRIw0lOR9HlolEEkVJFHEL14lEWMcBWzrfLnLIa-fL2Gyun2kj0C2J8uk2m_bXr0_D5jsiCrjxaiG5_DojxujNx75YVTjUeypD-xp0UATFQqc6VNjzckSGxZsA178jaS8-g6OZXglcWCl9EzS3AMxN6f39CbVsdHgKGgTUHSozqJLwyBCKLJfcznPEtxwVsUcn4KvdOtK6sJcmC4ibWF0S2O3gRVMXo9dUeLqFBzUtQMjcNH_EeSfxEGFiuR30C1ZM70v9T6GUD1N0PHhug__2F-DsN57m2ozwjTihQ../download [following]
    --2019-12-12 16:45:22--  https://public.boxcloud.com/d/1/b1!jD0baoEwhAs7CayO95UyhSNszeFfA_EqNIM7jM-ZixrgAG6UaRFMLlmYPaJ8YVAXdFfxrrHnQc0CJHNUY_6RT3QrlkhFFuq92WCeuhmH1Ltvb75ni-GO-K2Ltbni0HMOij3IAWXLEz4IqZLNeMjY1Ln2JyW8cLN_63RoCKDAkqv8okgYR_sB92M6_v2_Cuj3UvfzH4u8H6DKzeFch8_SfF7gubKOAFFmgzG1X4z6tXQBAgf6ezIfkS-K20C3Uw5PK-wxFvUVLIHi_xm9W4FaubPBkxbVlBUqB0WFhMLOlKty3KcLKptS67fy8xMpYQGVmhAUCUzkmOK33rKOBSdePbM9XeicLsreHA1xq41eYCmE_S_3-7YOMFMoS-3wjEjsatR-yL_wlHR7Xop2I4Qr9-1JEpBye_C8TtUjLeC7OPzvxvnZJ4r4nWtY1qH9dTJODRf2B9twvXLoE8HE-MR3noULJvEArKEn5d2dIDk1wSpNnHPkKhOhNTp1RP8fjTQXBtUs8vbymUAhDX2lhq0pR4mjbR9Rwphvvom0hpeeuBFEI1R-Wn4M7GnAzdf0cqNwTf2OwuMIq5spmd2DAVRe7QHSSiHEkpd4BDlXg2JTuKIhX4l4c1naAM_3g1zD9EdBWN3pdHamjn3G5Gtk3KuIoMKuOy_fnvIo-nKnrO3DE4J0bgim-GdTQHRbsxBGu2MRU9AVVfgw_yaf2c6Hd_NkAd0GFFJgyEelWVfLYu214F0njqpRY6n3hqDRC_XRSRPbEH7-Xv0hEBYZzZqH31DKc3jb2sTPsG-jaF4Sv7WR8qOtjD2Xk6uv-1MWXefO-pw9GsLmgNGhWX85gmEbAHIxqKA_oJq9pSi3bFNbpWrimvxFGIEcr42b6Utg1oUP35pL5_T8_04SkIhxDihD9gxe6EoYmhpM2kzbcUyViEnG_RUBIFCj7kDaQb71ODQFYPxTXMyb3c_PipXe_GGm_JtZCp_xJeOOZZAs4QGbRXNaznpTzzyA3r66sSycCpbl18rkQQybOoHo8XNhkdzCt-cB3Dw_dkN_YJQ8_BkUIiiij3yETRPEO-NTvHNX_pqWWMRgxNb__RGzoMiANsGyqFlFZmQZnS03gEx6zwa7LhPsTxI6797c1ichD8uXRIw0lOR9HlolEEkVJFHEL14lEWMcBWzrfLnLIa-fL2Gyun2kj0C2J8uk2m_bXr0_D5jsiCrjxaiG5_DojxujNx75YVTjUeypD-xp0UATFQqc6VNjzckSGxZsA178jaS8-g6OZXglcWCl9EzS3AMxN6f39CbVsdHgKGgTUHSozqJLwyBCKLJfcznPEtxwVsUcn4KvdOtK6sJcmC4ibWF0S2O3gRVMXo9dUeLqFBzUtQMjcNH_EeSfxEGFiuR30C1ZM70v9T6GUD1N0PHhug__2F-DsN57m2ozwjTihQ../download
    Resolving public.boxcloud.com (public.boxcloud.com)... 103.116.4.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|103.116.4.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 4547896018 (4.2G) [application/octet-stream]
    Saving to: ‘SRR6470906_S1_L001_R1_001.fastq.gz’
    
    SRR6470906_S1_L001_ 100%[===================>]   4.24G  15.9MB/s    in 4m 38s  
    
    2019-12-12 16:50:01 (15.6 MB/s) - ‘SRR6470906_S1_L001_R1_001.fastq.gz’ saved [4547896018/4547896018]
    
    --2019-12-12 16:50:03--  https://caltech.box.com/shared/static/63fh2xa5t82x7s74rqa0e2u2ur59y5ox.gz
    Resolving caltech.box.com (caltech.box.com)... 103.116.4.197
    Connecting to caltech.box.com (caltech.box.com)|103.116.4.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/63fh2xa5t82x7s74rqa0e2u2ur59y5ox.gz [following]
    --2019-12-12 16:50:03--  https://caltech.box.com/public/static/63fh2xa5t82x7s74rqa0e2u2ur59y5ox.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/63fh2xa5t82x7s74rqa0e2u2ur59y5ox.gz [following]
    --2019-12-12 16:50:04--  https://caltech.app.box.com/public/static/63fh2xa5t82x7s74rqa0e2u2ur59y5ox.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 103.116.4.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|103.116.4.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!IQRKpPtIYqhfyougjfYK79FUEIhKonLIzymhygPMsLCU9NejeJmaqCEIpnArmK_HhcOkAKi9KfMn8QjoRS1oz9JgPvc4dIP1jmr6XH6yHNoPtUOu-zH7QjxeAfr_aSbUESOQMljUQ-QwuQtkB-2ZuE24H58iNYURr0L1OAB6Wd_7DPeoqPQcuY9xlI4YESc8jmZZFcVKen7vSyx5p9vJY2SuzWrbt8JTiJ5By7APStib5HVw6vFbitaxU-ibdN4uowrGgfzDRSLmyIA6CwaePNehpPCHjfQfF0rieN19M3pxAUG8vEgsvG2VrMUCJ1aZMs1sTr9c6y9m1AOIDAUJheqfqLvfhbs5tsGskwVQvEbOn3cj8oTz2sxSOr01axHwofC7jepknJHS-quJVGFUgp0DfAYIbvTYQ5uQfC0iXxcn4bnPLkYBRHh29tPSPGN2unBDD2TAludj6xhXeeX47LNGa62JkuUru6DdlpLdZuJPpnt3VoowVE8OBpStlyt2p7rEdZUfPPAt9_AW16DZYScA7e3DR7daODJrN9n85lQTKT279M2XeeomRn7s84ZfqxUXshh9kHBx06NqmTM-IOlsFfN7epF27AwfQGAu1tAwCkNrVDhEKWZnllMCjUAMhamqDn_CkGdlrlSIT0gYeGFwttAywRbXbXv-FOlOx6sWhnUWSSxQy6FeG8dgsbfGyLeB0MAC-I2o_sj_2bLgD0gb3b9zatzPcuIm4CjvYPzEdYbI5m-AfGMfhRfTWJchurBxbadhrEG77NRjILyQsmcm4J7VlrUrf1nv8vUADLRbHYsclkbYKKuvLz6o34J9bspPNatCrum06h750_YyDsCEC2D27iap1k9r_aPkzAba1yDlWAu-oq6yPMN5ROtoCqnGTJD2Pl7xzrxFUjMaVF31lHqodRQXHem9VctLmeZthACgi8GNZJr5Jvil9zXtCTHX39xcJF0uBOvrC7O5-4dDj5qB0jSl0DsAWPLhHG1v0s1YdzuySAvt-DSyfrqeuRD9S9V1eUxiC_wiP2IUNGJIi3Ve3aiHA4t8D4iC4eSzkkwKSAEtL60I6KNg4lZdrFd-i2oVuvkfZC0IlV7zVaCP3cgIu851WRyDx7gYy4CXcOdjvfgbHmKXoHNroLF1dIbV5Edht3MN5cfRO7T8TVBPv15Fq4TPxBnJa5zEtIF_UoyqczbQvFQH0sTbDyAHeW-N0Y5-LkA-jAdO0oE-B3vaPq7KIZK-6w09pFXiQ_5Sq-9WUy_yGEfOGCDoRwybfzbtWniQcVTpoQatU1GGlyk2EOmPF8iU6WDdDIPL5WpOdba8w-48w6_L2D-SSNYFFmLkN6LhvZ12_1Dt__azAM85EBahsFiVUIWSrn7OIPwWBiuLXc_XqEB0e6IbF7FZM57YeArhrukqmuK_yk5u2vKC8S9ofQ../download [following]
    --2019-12-12 16:50:04--  https://public.boxcloud.com/d/1/b1!IQRKpPtIYqhfyougjfYK79FUEIhKonLIzymhygPMsLCU9NejeJmaqCEIpnArmK_HhcOkAKi9KfMn8QjoRS1oz9JgPvc4dIP1jmr6XH6yHNoPtUOu-zH7QjxeAfr_aSbUESOQMljUQ-QwuQtkB-2ZuE24H58iNYURr0L1OAB6Wd_7DPeoqPQcuY9xlI4YESc8jmZZFcVKen7vSyx5p9vJY2SuzWrbt8JTiJ5By7APStib5HVw6vFbitaxU-ibdN4uowrGgfzDRSLmyIA6CwaePNehpPCHjfQfF0rieN19M3pxAUG8vEgsvG2VrMUCJ1aZMs1sTr9c6y9m1AOIDAUJheqfqLvfhbs5tsGskwVQvEbOn3cj8oTz2sxSOr01axHwofC7jepknJHS-quJVGFUgp0DfAYIbvTYQ5uQfC0iXxcn4bnPLkYBRHh29tPSPGN2unBDD2TAludj6xhXeeX47LNGa62JkuUru6DdlpLdZuJPpnt3VoowVE8OBpStlyt2p7rEdZUfPPAt9_AW16DZYScA7e3DR7daODJrN9n85lQTKT279M2XeeomRn7s84ZfqxUXshh9kHBx06NqmTM-IOlsFfN7epF27AwfQGAu1tAwCkNrVDhEKWZnllMCjUAMhamqDn_CkGdlrlSIT0gYeGFwttAywRbXbXv-FOlOx6sWhnUWSSxQy6FeG8dgsbfGyLeB0MAC-I2o_sj_2bLgD0gb3b9zatzPcuIm4CjvYPzEdYbI5m-AfGMfhRfTWJchurBxbadhrEG77NRjILyQsmcm4J7VlrUrf1nv8vUADLRbHYsclkbYKKuvLz6o34J9bspPNatCrum06h750_YyDsCEC2D27iap1k9r_aPkzAba1yDlWAu-oq6yPMN5ROtoCqnGTJD2Pl7xzrxFUjMaVF31lHqodRQXHem9VctLmeZthACgi8GNZJr5Jvil9zXtCTHX39xcJF0uBOvrC7O5-4dDj5qB0jSl0DsAWPLhHG1v0s1YdzuySAvt-DSyfrqeuRD9S9V1eUxiC_wiP2IUNGJIi3Ve3aiHA4t8D4iC4eSzkkwKSAEtL60I6KNg4lZdrFd-i2oVuvkfZC0IlV7zVaCP3cgIu851WRyDx7gYy4CXcOdjvfgbHmKXoHNroLF1dIbV5Edht3MN5cfRO7T8TVBPv15Fq4TPxBnJa5zEtIF_UoyqczbQvFQH0sTbDyAHeW-N0Y5-LkA-jAdO0oE-B3vaPq7KIZK-6w09pFXiQ_5Sq-9WUy_yGEfOGCDoRwybfzbtWniQcVTpoQatU1GGlyk2EOmPF8iU6WDdDIPL5WpOdba8w-48w6_L2D-SSNYFFmLkN6LhvZ12_1Dt__azAM85EBahsFiVUIWSrn7OIPwWBiuLXc_XqEB0e6IbF7FZM57YeArhrukqmuK_yk5u2vKC8S9ofQ../download
    Resolving public.boxcloud.com (public.boxcloud.com)... 103.116.4.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|103.116.4.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 8051434084 (7.5G) [application/octet-stream]
    Saving to: ‘SRR6470906_S1_L001_R2_001.fastq.gz’
    
    SRR6470906_S1_L001_ 100%[===================>]   7.50G  13.8MB/s    in 9m 14s  
    
    2019-12-12 16:59:19 (13.9 MB/s) - ‘SRR6470906_S1_L001_R2_001.fastq.gz’ saved [8051434084/8051434084]
    
    --2019-12-12 16:59:21--  https://caltech.box.com/shared/static/zqi3durukillaw1pbns1kd1lowyfg5qk.gz
    Resolving caltech.box.com (caltech.box.com)... 103.116.4.197
    Connecting to caltech.box.com (caltech.box.com)|103.116.4.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/zqi3durukillaw1pbns1kd1lowyfg5qk.gz [following]
    --2019-12-12 16:59:22--  https://caltech.box.com/public/static/zqi3durukillaw1pbns1kd1lowyfg5qk.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/zqi3durukillaw1pbns1kd1lowyfg5qk.gz [following]
    --2019-12-12 16:59:22--  https://caltech.app.box.com/public/static/zqi3durukillaw1pbns1kd1lowyfg5qk.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 103.116.4.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|103.116.4.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!GOJlDv0dHx7IF1T0A37vi_RMeDchTdIRIbmHryH43jHU4uSYmqsYBNq2UVYgHo4nHavUlf9BczDrOJo9v7aZ2ApTifBKTx0xE9QVSAf2jAhqtbFGufK9QhmNpu1yVhsEZdF1WqrxgI_FkGwLqQ3Bx8ZuFH2o4DAQc4MLS0HnUhHDiCGUoV6Z8C-kAb0f_2YEcQgSfjaVsBNTXBzrwnEYPjzqDs0eesG3fQN_99O-tTt6MKRSWZHkmHiCSEA8nToRCPoF3c72gmxAduzKpZvudzj_avmyH5EwdYjOOy63k61fr3eJVQdDivNaTm88_9NKbzj2-HvlX9MV2a7FmFFAmkKtvZtoc-JdzIXXGb241Vs5ewwhOIQX3pikJxvHC_Wo1SUsWHyf0m35348Z0ef3VS_9rHDVGvTdtG-waRcSOfAqpWmxaRKs_Qg7IYVYEvAfJP3C6hFT3QglnrTSkXc2FHwy4uDq76RonZNi033xdhlPIbFUTXpEieUVUYjpS6YiwfxvSX1K9C3L5J71C02Vtt_3M3It4uTeKdzxzpsUhKdIMyAKRhE1ZLiLbEStihDg33gGEzz65ZP6tWUQPwUYomj2N9HL_wZ8-9bD7CxDVlA8sQ6yc0HCO2P9UUTfaJuY_q9lFo2Iz-e4MxOLBT9SnRirzGAqmxi0lq4xArN5euZ4ltr-yA-VXiR3-ai3hdxXl23sYz-5bz6k5629zb87fCDwb4MMIuMV0pw00_N7RvS6G9t6knEAZAiJZshI5p21Z0nogmPQcj9LodjnhGZXkpaUS_J0eP3LxK8UeUs55y9H07rYXdQlBZL27XBNP2vRhCChImZsoEzQAw5Yqq4vpR064mSj00Luulo-tF62ihTToc4icGQmqNPfagFkOxkq6YvjudgQybkL7dI-l7nXzcQ2rld_jU6svLY5e3v-kB_oRMIxdPwCd3lwNr6D7I0M21jTngGVFhkQb132dsjb9zSge50uekCYSOKb77Rw1JztzQxurG3LH-t85dnfA-26xeHV-EL7Oi3350GwNlqV74yGN540ZeLKYaFCdieR_tpoJFz0NyxFsaAQE2I2BMoyeVqKIiEP4_eCn0U2ihA7WG01bfXMCYHhtLw4DAoOtCJzUfweCQw9m2VZLaiMXtY9ymTFVthvahyBWVH4njLCNVKwLZEk52aqdsXsv7yU2KpaE4Mg7MtEJPDatPjxIHoZLXEuNSSRMygmrC0GTsSAvKfk1zJVZYh_tHsHk_yBaXi6TIZDj9b4CcQSwiF3j6kS7--6WDivkAFx_Sut55iOQoEWoJAvrChKSOtp2F3d0_YBbKlc9M_5Y1L4Yg9tGRVtTop6cCyysoNbJB5EKVKX9MtcxJdfl9aA_14uz_-59i22FiOYpHTFAUw4iYWNmpk5sWu6Uu_HoKag498Y5xaMvmG5fm9n-g../download [following]
    --2019-12-12 16:59:23--  https://public.boxcloud.com/d/1/b1!GOJlDv0dHx7IF1T0A37vi_RMeDchTdIRIbmHryH43jHU4uSYmqsYBNq2UVYgHo4nHavUlf9BczDrOJo9v7aZ2ApTifBKTx0xE9QVSAf2jAhqtbFGufK9QhmNpu1yVhsEZdF1WqrxgI_FkGwLqQ3Bx8ZuFH2o4DAQc4MLS0HnUhHDiCGUoV6Z8C-kAb0f_2YEcQgSfjaVsBNTXBzrwnEYPjzqDs0eesG3fQN_99O-tTt6MKRSWZHkmHiCSEA8nToRCPoF3c72gmxAduzKpZvudzj_avmyH5EwdYjOOy63k61fr3eJVQdDivNaTm88_9NKbzj2-HvlX9MV2a7FmFFAmkKtvZtoc-JdzIXXGb241Vs5ewwhOIQX3pikJxvHC_Wo1SUsWHyf0m35348Z0ef3VS_9rHDVGvTdtG-waRcSOfAqpWmxaRKs_Qg7IYVYEvAfJP3C6hFT3QglnrTSkXc2FHwy4uDq76RonZNi033xdhlPIbFUTXpEieUVUYjpS6YiwfxvSX1K9C3L5J71C02Vtt_3M3It4uTeKdzxzpsUhKdIMyAKRhE1ZLiLbEStihDg33gGEzz65ZP6tWUQPwUYomj2N9HL_wZ8-9bD7CxDVlA8sQ6yc0HCO2P9UUTfaJuY_q9lFo2Iz-e4MxOLBT9SnRirzGAqmxi0lq4xArN5euZ4ltr-yA-VXiR3-ai3hdxXl23sYz-5bz6k5629zb87fCDwb4MMIuMV0pw00_N7RvS6G9t6knEAZAiJZshI5p21Z0nogmPQcj9LodjnhGZXkpaUS_J0eP3LxK8UeUs55y9H07rYXdQlBZL27XBNP2vRhCChImZsoEzQAw5Yqq4vpR064mSj00Luulo-tF62ihTToc4icGQmqNPfagFkOxkq6YvjudgQybkL7dI-l7nXzcQ2rld_jU6svLY5e3v-kB_oRMIxdPwCd3lwNr6D7I0M21jTngGVFhkQb132dsjb9zSge50uekCYSOKb77Rw1JztzQxurG3LH-t85dnfA-26xeHV-EL7Oi3350GwNlqV74yGN540ZeLKYaFCdieR_tpoJFz0NyxFsaAQE2I2BMoyeVqKIiEP4_eCn0U2ihA7WG01bfXMCYHhtLw4DAoOtCJzUfweCQw9m2VZLaiMXtY9ymTFVthvahyBWVH4njLCNVKwLZEk52aqdsXsv7yU2KpaE4Mg7MtEJPDatPjxIHoZLXEuNSSRMygmrC0GTsSAvKfk1zJVZYh_tHsHk_yBaXi6TIZDj9b4CcQSwiF3j6kS7--6WDivkAFx_Sut55iOQoEWoJAvrChKSOtp2F3d0_YBbKlc9M_5Y1L4Yg9tGRVtTop6cCyysoNbJB5EKVKX9MtcxJdfl9aA_14uz_-59i22FiOYpHTFAUw4iYWNmpk5sWu6Uu_HoKag498Y5xaMvmG5fm9n-g../download
    Resolving public.boxcloud.com (public.boxcloud.com)... 103.116.4.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|103.116.4.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 3473816553 (3.2G) [application/octet-stream]
    Saving to: ‘SRR6470906_S1_L002_R1_001.fastq.gz’
    
    SRR6470906_S1_L002_ 100%[===================>]   3.23G  14.8MB/s    in 3m 50s  
    
    2019-12-12 17:03:13 (14.4 MB/s) - ‘SRR6470906_S1_L002_R1_001.fastq.gz’ saved [3473816553/3473816553]
    
    --2019-12-12 17:03:16--  https://caltech.box.com/shared/static/i56qojfz41ns1kw9z86sla0vawsch96t.gz
    Resolving caltech.box.com (caltech.box.com)... 103.116.4.197
    Connecting to caltech.box.com (caltech.box.com)|103.116.4.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/i56qojfz41ns1kw9z86sla0vawsch96t.gz [following]
    --2019-12-12 17:03:17--  https://caltech.box.com/public/static/i56qojfz41ns1kw9z86sla0vawsch96t.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/i56qojfz41ns1kw9z86sla0vawsch96t.gz [following]
    --2019-12-12 17:03:17--  https://caltech.app.box.com/public/static/i56qojfz41ns1kw9z86sla0vawsch96t.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 103.116.4.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|103.116.4.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!R6y7a5bIJuqYMSWPoHnEy-l6dd72cFwhiG5GMUeOmCo3_B5BccUZA-SpS3dRR8Zekhj0VXOjdARtBlPqDNC8Y4ol4fYDVjJhh9XFkGLQ-1Jzjr3kSxEWbjPkbI_6gDfAmHJvxZifs5Dhuy7svYvomAgSezzR6cFIBwQounJ6mHB5S68EQy6Mh4xs_aTunCVHqDjaTn1GPZfoeDQ297sblMskuI2Aq3mdVCwkrg5suolP7FHHdz7YwoARSDx3mE3o1IdwntBZGueppadjfdydZR_icMVAywN3zE_AFCNKRAPEVqbEpEAe9uQrrO9VtmgC9hDEo6wzObx1yIpQ3YLerB6i8tVXubO9mWHj6D-ofNRaDNr0TCMkwBxtCD8smWDiDw3avJOKy1YcqOco5PgP0YnjYFuuPxQlNRD-jTuXZNmallPtZDteXuTv0kFCV2K0iWZwHgpyqZBmLoQk5OiFBjLOZ28pz0YenhK5k60bRxKksr6MD6PnTWNzpgIXGEfPrrtz0V2D7i_b-4xvNep23Ws0ExS6Yrb0gYDDNi8ewWvDgunmmvZxvW8UfqAKSFg4FNNu3_9wSgsutfidU-2Q5Qgd3Ozx9WU7czC8Cf_PJWlD6pJKnXVKGVrYulXQs2l6bYFOncJWYxGJiojwaL-rFDC45EupjnosxnacedbsdT1qmO3cY8hwYvceRzZ9_HZErKGPnELgPgUHijFcd2z7rvxmpztQePEocCwaucO3byUJ_zkjfeNWnsmHLn8s1oFF6gjA1-cQKSOggR7hcqWwsFgtFdBbeOy9TAAnYpiMhtH2s_KnoRiA_aJDOJRT0huK3gB-YT-4U0EsBGLkYw6np7txv6OXT7TDti98EE4XvQGFeu9-pSGFSvUirbt02uL0RJtQX77TEZH9-_18bfH8YlatHgfuc1fqo-pDoOkccIBekbTTCz5LMT-zhBnb8_GgyVTJlBsxpeP46QgiaAqMyNCHC_sXzdYBPeaQQtVOOSp1TSdGIbcqKxzgn5ToELW6q623X5BgT1-P-0E2Rqyn90o5EyCjtaiUOW6maWcr-jGiydXWMcypp4sE-5OyT8hWT4Ipu1UVVFi0U1QtAfKOZ_zQD5P1OY3RRTQdqVUGUksu42t4KPJneopa73QgMXTsSdXl2RnSG2yFuqbHxFQT3SwOrz1gPaq9JPGv0PFqMU-7QNThEwiWTo7UTNrR-FG-xEFPC8jbqNaSX7lkrstLXkFbGjy4sTHqgpn2RSVBfc340wKvMzRraWtoQSer990foD--oyoCdm7bz8dWbKDCSN7_zmiu0CQShK6M4EWstjyfBZdN8F3BKojjcIba9EwBvclmH_u_MBb6TRbDfTnbCbljYi43vXzFwb_uPYLCdiVO6lbr3t6R1u_LEo5OuEAqWXAkTbqgEzBs27uhPS-IWvTruGP1FQ../download [following]
    --2019-12-12 17:03:18--  https://public.boxcloud.com/d/1/b1!R6y7a5bIJuqYMSWPoHnEy-l6dd72cFwhiG5GMUeOmCo3_B5BccUZA-SpS3dRR8Zekhj0VXOjdARtBlPqDNC8Y4ol4fYDVjJhh9XFkGLQ-1Jzjr3kSxEWbjPkbI_6gDfAmHJvxZifs5Dhuy7svYvomAgSezzR6cFIBwQounJ6mHB5S68EQy6Mh4xs_aTunCVHqDjaTn1GPZfoeDQ297sblMskuI2Aq3mdVCwkrg5suolP7FHHdz7YwoARSDx3mE3o1IdwntBZGueppadjfdydZR_icMVAywN3zE_AFCNKRAPEVqbEpEAe9uQrrO9VtmgC9hDEo6wzObx1yIpQ3YLerB6i8tVXubO9mWHj6D-ofNRaDNr0TCMkwBxtCD8smWDiDw3avJOKy1YcqOco5PgP0YnjYFuuPxQlNRD-jTuXZNmallPtZDteXuTv0kFCV2K0iWZwHgpyqZBmLoQk5OiFBjLOZ28pz0YenhK5k60bRxKksr6MD6PnTWNzpgIXGEfPrrtz0V2D7i_b-4xvNep23Ws0ExS6Yrb0gYDDNi8ewWvDgunmmvZxvW8UfqAKSFg4FNNu3_9wSgsutfidU-2Q5Qgd3Ozx9WU7czC8Cf_PJWlD6pJKnXVKGVrYulXQs2l6bYFOncJWYxGJiojwaL-rFDC45EupjnosxnacedbsdT1qmO3cY8hwYvceRzZ9_HZErKGPnELgPgUHijFcd2z7rvxmpztQePEocCwaucO3byUJ_zkjfeNWnsmHLn8s1oFF6gjA1-cQKSOggR7hcqWwsFgtFdBbeOy9TAAnYpiMhtH2s_KnoRiA_aJDOJRT0huK3gB-YT-4U0EsBGLkYw6np7txv6OXT7TDti98EE4XvQGFeu9-pSGFSvUirbt02uL0RJtQX77TEZH9-_18bfH8YlatHgfuc1fqo-pDoOkccIBekbTTCz5LMT-zhBnb8_GgyVTJlBsxpeP46QgiaAqMyNCHC_sXzdYBPeaQQtVOOSp1TSdGIbcqKxzgn5ToELW6q623X5BgT1-P-0E2Rqyn90o5EyCjtaiUOW6maWcr-jGiydXWMcypp4sE-5OyT8hWT4Ipu1UVVFi0U1QtAfKOZ_zQD5P1OY3RRTQdqVUGUksu42t4KPJneopa73QgMXTsSdXl2RnSG2yFuqbHxFQT3SwOrz1gPaq9JPGv0PFqMU-7QNThEwiWTo7UTNrR-FG-xEFPC8jbqNaSX7lkrstLXkFbGjy4sTHqgpn2RSVBfc340wKvMzRraWtoQSer990foD--oyoCdm7bz8dWbKDCSN7_zmiu0CQShK6M4EWstjyfBZdN8F3BKojjcIba9EwBvclmH_u_MBb6TRbDfTnbCbljYi43vXzFwb_uPYLCdiVO6lbr3t6R1u_LEo5OuEAqWXAkTbqgEzBs27uhPS-IWvTruGP1FQ../download
    Resolving public.boxcloud.com (public.boxcloud.com)... 103.116.4.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|103.116.4.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 6150590468 (5.7G) [application/octet-stream]
    Saving to: ‘SRR6470906_S1_L002_R2_001.fastq.gz’
    
    SRR6470906_S1_L002_ 100%[===================>]   5.73G  16.0MB/s    in 6m 17s  
    
    2019-12-12 17:09:36 (15.6 MB/s) - ‘SRR6470906_S1_L002_R2_001.fastq.gz’ saved [6150590468/6150590468]
    
    --2019-12-12 17:09:38--  https://caltech.box.com/shared/static/vrditbbk38tw3f61fwpg504vcc5x09ci.gz
    Resolving caltech.box.com (caltech.box.com)... 103.116.4.197
    Connecting to caltech.box.com (caltech.box.com)|103.116.4.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/vrditbbk38tw3f61fwpg504vcc5x09ci.gz [following]
    --2019-12-12 17:09:38--  https://caltech.box.com/public/static/vrditbbk38tw3f61fwpg504vcc5x09ci.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/vrditbbk38tw3f61fwpg504vcc5x09ci.gz [following]
    --2019-12-12 17:09:39--  https://caltech.app.box.com/public/static/vrditbbk38tw3f61fwpg504vcc5x09ci.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 103.116.4.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|103.116.4.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!mEw7jD3SYx2JkfpPGUANWAkasmzfYG_4IsWs_1IIcIove2w1TUwrVVpSqlNM8xPs-NiiaFcTJOO191k11vZNt5UaleiAtPOJbdCvrwABOUgpE17vVBhKIsjOpp_DmmYeivZaVBJMh_GQGiVKzkhVwTkbzm5s-EZunbREKZsvTfjrzWTUv_xaHoT5deMNloKiCw--jQoahlEPn4Ng3uOSpHTlGfjcx-LcwzZS9mxrkfsW7PJGpaTn3PiDw-u6PuI4RuN9BmviJz0fqJs2vLOdKuGHOUNnV3p3DgnHX-NxFtYxBXf2m_xC-iYN3hXPbQF1UvruvypWCzIKLMzkyJRGTlgMENXff9YtLB04B4bATBJw9xO8YLEA1UaQ-2DYdi4pRY5tSM_u_KwfYhTsBj0PzUOjSUe3plTIOPIeriZAB02SdU2mGMuQOQrzNy0FXgFYI3nMmUPdzTGYA4x5vt6dP99T9wwyRromJMLtfY-sjiK3j_WgUjxvP6Wpm7IM8vGzls5ZrJ85LktkQDyU7sOZLojFgjITbaPSNbJDtcaiI23iUNeJELerAM-opkV8CbgV4ohHWrgQIUyKm45mu0K8kbGjo6Ve40mgfDb1CdsX7zBWPEDfd9wwXLq3c3C1DUSuW0_9kvZSreQ6F-lNaGwLSD62Y4ruh7mRwZMHr90Y3nY0JnsXbaZMDXCpDLRuCMu_oTOR2LRY3fU2Ez6s0m4gtlgRFzsvwtqXAHCHO8Rj4kXq0HrpYljhxiF1c8q-PZqvEPcf5DmnT4zpEr9KG5gpqRlCwwfJWl1dJ-3U4SxvXATHQtwqM1LLLO2hDv8GmuQ05tERAlWHRrx0zOT9cmjSAs63b-6BfApttPjr55oxKp59ouozkWBbNllndr_qnML8faaoak3Yy4wq8NfWJeRIejfTEaNg_Gv5trktY99C4sPEcxHjeRcY6HWhR0Z1YhpIhcynayVMEAhgUGNRbtEhLF-hx6ZqJUKv8gBKmxmKwZoNjmoEJiOTDULD9tl3cp8Tlj9FjraA7Fl3BCVWo4o56fCCgqs8h9txsYaYbHcyLMtA9LdWwpgLf9EcOEObmSMsi9Dt-ORHTpkMa7jdiBOF1uuczrmK58fj_Mi46IdNjiLB3_595-BK2ClzeWsu7Os-RznkS-BNc1gWFgB5T_1wsxOTpJ3X0RKG3LHRQHgd2yGpXqStAS_j7Z1ze0ooWMO1vFbwZdca05sq-DnIvNy6mtY4Jt4DKconl562MC5FEI4QKtW_FSZPFbWjC7Nioo9CTYJbRbPM6ADHvpjL4RrW2T0MumHF9V1TgF7wDLo79i9N-3BBxvWFruNTi-iMJ30zvvQ0k2Xd1wrbeSflfhS4mtZbvf6nOJZr5COaSsqcrT8BzViSwr09bocU54hRTB_5MFg0NilGApeu3eL20kBxNsf_zjehrVU./download [following]
    --2019-12-12 17:09:39--  https://public.boxcloud.com/d/1/b1!mEw7jD3SYx2JkfpPGUANWAkasmzfYG_4IsWs_1IIcIove2w1TUwrVVpSqlNM8xPs-NiiaFcTJOO191k11vZNt5UaleiAtPOJbdCvrwABOUgpE17vVBhKIsjOpp_DmmYeivZaVBJMh_GQGiVKzkhVwTkbzm5s-EZunbREKZsvTfjrzWTUv_xaHoT5deMNloKiCw--jQoahlEPn4Ng3uOSpHTlGfjcx-LcwzZS9mxrkfsW7PJGpaTn3PiDw-u6PuI4RuN9BmviJz0fqJs2vLOdKuGHOUNnV3p3DgnHX-NxFtYxBXf2m_xC-iYN3hXPbQF1UvruvypWCzIKLMzkyJRGTlgMENXff9YtLB04B4bATBJw9xO8YLEA1UaQ-2DYdi4pRY5tSM_u_KwfYhTsBj0PzUOjSUe3plTIOPIeriZAB02SdU2mGMuQOQrzNy0FXgFYI3nMmUPdzTGYA4x5vt6dP99T9wwyRromJMLtfY-sjiK3j_WgUjxvP6Wpm7IM8vGzls5ZrJ85LktkQDyU7sOZLojFgjITbaPSNbJDtcaiI23iUNeJELerAM-opkV8CbgV4ohHWrgQIUyKm45mu0K8kbGjo6Ve40mgfDb1CdsX7zBWPEDfd9wwXLq3c3C1DUSuW0_9kvZSreQ6F-lNaGwLSD62Y4ruh7mRwZMHr90Y3nY0JnsXbaZMDXCpDLRuCMu_oTOR2LRY3fU2Ez6s0m4gtlgRFzsvwtqXAHCHO8Rj4kXq0HrpYljhxiF1c8q-PZqvEPcf5DmnT4zpEr9KG5gpqRlCwwfJWl1dJ-3U4SxvXATHQtwqM1LLLO2hDv8GmuQ05tERAlWHRrx0zOT9cmjSAs63b-6BfApttPjr55oxKp59ouozkWBbNllndr_qnML8faaoak3Yy4wq8NfWJeRIejfTEaNg_Gv5trktY99C4sPEcxHjeRcY6HWhR0Z1YhpIhcynayVMEAhgUGNRbtEhLF-hx6ZqJUKv8gBKmxmKwZoNjmoEJiOTDULD9tl3cp8Tlj9FjraA7Fl3BCVWo4o56fCCgqs8h9txsYaYbHcyLMtA9LdWwpgLf9EcOEObmSMsi9Dt-ORHTpkMa7jdiBOF1uuczrmK58fj_Mi46IdNjiLB3_595-BK2ClzeWsu7Os-RznkS-BNc1gWFgB5T_1wsxOTpJ3X0RKG3LHRQHgd2yGpXqStAS_j7Z1ze0ooWMO1vFbwZdca05sq-DnIvNy6mtY4Jt4DKconl562MC5FEI4QKtW_FSZPFbWjC7Nioo9CTYJbRbPM6ADHvpjL4RrW2T0MumHF9V1TgF7wDLo79i9N-3BBxvWFruNTi-iMJ30zvvQ0k2Xd1wrbeSflfhS4mtZbvf6nOJZr5COaSsqcrT8BzViSwr09bocU54hRTB_5MFg0NilGApeu3eL20kBxNsf_zjehrVU./download
    Resolving public.boxcloud.com (public.boxcloud.com)... 103.116.4.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|103.116.4.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 3828633558 (3.6G) [application/octet-stream]
    Saving to: ‘SRR6470907_S1_L001_R1_001.fastq.gz’
    
    SRR6470907_S1_L001_ 100%[===================>]   3.57G  15.7MB/s    in 3m 58s  
    
    2019-12-12 17:13:38 (15.4 MB/s) - ‘SRR6470907_S1_L001_R1_001.fastq.gz’ saved [3828633558/3828633558]
    
    --2019-12-12 17:13:40--  https://caltech.box.com/shared/static/8ud3otwztjeqlmjctbu1fw7hg3k56ejr.gz
    Resolving caltech.box.com (caltech.box.com)... 103.116.4.197
    Connecting to caltech.box.com (caltech.box.com)|103.116.4.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/8ud3otwztjeqlmjctbu1fw7hg3k56ejr.gz [following]
    --2019-12-12 17:13:40--  https://caltech.box.com/public/static/8ud3otwztjeqlmjctbu1fw7hg3k56ejr.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/8ud3otwztjeqlmjctbu1fw7hg3k56ejr.gz [following]
    --2019-12-12 17:13:40--  https://caltech.app.box.com/public/static/8ud3otwztjeqlmjctbu1fw7hg3k56ejr.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 103.116.4.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|103.116.4.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!woL8LjzIjNqeJ9OSHk7-Xn9zi3K2VdvRVvUbaVnaQJf_poiZY3cA1FkCz-mNcDiv8inLCBr9HxopHLD7qtCpM0ZI3CWVA3sYKKYP146oeDkFtwrnKSX53aUm65dPCv2m1Q3IEkxI5SwD10xXsCWO1TJ91FDo84VieVTyqg5LbREx06zlJcYWvhwayUrSoFehJj1rSSmB27EZC9CgsphjiYUnoyiNyXegv0hy63_3Kwr-G1SKgbA4AO5L5nYvv_murpqwxalzTQgtHpmzTGcCVVERxQA4FULICpd-VmCxERdwUjhvCr5Uiu6DV499heqeH7IEPs-D2VH5VrdQY8gjgQ9_nK3T6Z1UcHPjgcewfqyJpvFVl1NQjTDmsBMUOBdftfIzhgRwJ_Iny6mrwTRgal1kYhwUMrbq6dgUy5-XhzodwXQrQItFDoYQeDKgLShOZMOs0P9hH_QJmEAEtuG5ax6KLfGpJzjoQJUs0SEhEzOeW9SltM3n0Mz09botMDcDYHWPYPybjAK95uFD0TdIefcMv_BVERwD2vG6GQLivK1mZ6ARg9xS6zbYF7KkffDE6kDUpeWkLkl_CKT0_0p9z5NMWPwGitM0NSOzvV8103q5LyUNboOBc18XKa-H8fXYFcAhZLeqsEkQ0rhbPKgYnB3lm98sSQBGFW8FlgFvTCGf2NYhazjsD4ZKi11duuGhHXrRYr7W1UBfdW7LTZObQqvnyvXQ0P9va2wAW6rkRbs8Oqbb8nlbnaJ9hQKzTRHR3hSl-6K98Ocw5k-B5d7RgzMRRLFsuRDEjZvSYf9EXY8bXzI_xjeSD7-Bs1vbSf95gC9fnxb2hhdQQe_fZYeWYP3OCqstzAI67RDD8dVOLsQ122CpmJRfQ7Uvl0jci0XFiJK--GqfjQ3LwboZfcfS30tNOsDmYFkmwi4wF3P9nrEFaN12Mh9J5IH2-W44MPZx_vpS-LRDqRAES3VYtOHwR245hp2Hp27rZebkeEeKXh4uL0C38AAd-FaSXoMswUPyuXj963O5HY3F8dFkaBU2LsHhG5LMLlOMxmwFSWgBOsrvmnuCmDxar2fHqh-VDTRYtUI3aEBYX9zv6JGizpzqP51VdgOg2_3v75pkLW0z69Lkmm4rQVig9rlrvz-5WMS5Ess4di2M_PynKqgJRNCIi6ZX3HZNY9w3M-0GTvASRwVzfGv_FbTCkXZhbWGCSfqjPOWEgcDsjP3sJAZloat5NhJjGKZHE8gXiDwER1rKmP_xPWylu-DvPzVbnQxAUHlQNRwEn1HeCc0A5oLYihZswtgxj63nnRCiQC-J71hnhrR6RH8naR3GfcxQ66mSHcM6FvXw9OLevcTikF1Ryf_4s-V2HL8ODiShjroat2OzvLjrYmzazUrzZa_Bw-TEDqP-_eQG20t1Y9rnMKG5l4ayIa6CbduANyM./download [following]
    --2019-12-12 17:13:41--  https://public.boxcloud.com/d/1/b1!woL8LjzIjNqeJ9OSHk7-Xn9zi3K2VdvRVvUbaVnaQJf_poiZY3cA1FkCz-mNcDiv8inLCBr9HxopHLD7qtCpM0ZI3CWVA3sYKKYP146oeDkFtwrnKSX53aUm65dPCv2m1Q3IEkxI5SwD10xXsCWO1TJ91FDo84VieVTyqg5LbREx06zlJcYWvhwayUrSoFehJj1rSSmB27EZC9CgsphjiYUnoyiNyXegv0hy63_3Kwr-G1SKgbA4AO5L5nYvv_murpqwxalzTQgtHpmzTGcCVVERxQA4FULICpd-VmCxERdwUjhvCr5Uiu6DV499heqeH7IEPs-D2VH5VrdQY8gjgQ9_nK3T6Z1UcHPjgcewfqyJpvFVl1NQjTDmsBMUOBdftfIzhgRwJ_Iny6mrwTRgal1kYhwUMrbq6dgUy5-XhzodwXQrQItFDoYQeDKgLShOZMOs0P9hH_QJmEAEtuG5ax6KLfGpJzjoQJUs0SEhEzOeW9SltM3n0Mz09botMDcDYHWPYPybjAK95uFD0TdIefcMv_BVERwD2vG6GQLivK1mZ6ARg9xS6zbYF7KkffDE6kDUpeWkLkl_CKT0_0p9z5NMWPwGitM0NSOzvV8103q5LyUNboOBc18XKa-H8fXYFcAhZLeqsEkQ0rhbPKgYnB3lm98sSQBGFW8FlgFvTCGf2NYhazjsD4ZKi11duuGhHXrRYr7W1UBfdW7LTZObQqvnyvXQ0P9va2wAW6rkRbs8Oqbb8nlbnaJ9hQKzTRHR3hSl-6K98Ocw5k-B5d7RgzMRRLFsuRDEjZvSYf9EXY8bXzI_xjeSD7-Bs1vbSf95gC9fnxb2hhdQQe_fZYeWYP3OCqstzAI67RDD8dVOLsQ122CpmJRfQ7Uvl0jci0XFiJK--GqfjQ3LwboZfcfS30tNOsDmYFkmwi4wF3P9nrEFaN12Mh9J5IH2-W44MPZx_vpS-LRDqRAES3VYtOHwR245hp2Hp27rZebkeEeKXh4uL0C38AAd-FaSXoMswUPyuXj963O5HY3F8dFkaBU2LsHhG5LMLlOMxmwFSWgBOsrvmnuCmDxar2fHqh-VDTRYtUI3aEBYX9zv6JGizpzqP51VdgOg2_3v75pkLW0z69Lkmm4rQVig9rlrvz-5WMS5Ess4di2M_PynKqgJRNCIi6ZX3HZNY9w3M-0GTvASRwVzfGv_FbTCkXZhbWGCSfqjPOWEgcDsjP3sJAZloat5NhJjGKZHE8gXiDwER1rKmP_xPWylu-DvPzVbnQxAUHlQNRwEn1HeCc0A5oLYihZswtgxj63nnRCiQC-J71hnhrR6RH8naR3GfcxQ66mSHcM6FvXw9OLevcTikF1Ryf_4s-V2HL8ODiShjroat2OzvLjrYmzazUrzZa_Bw-TEDqP-_eQG20t1Y9rnMKG5l4ayIa6CbduANyM./download
    Resolving public.boxcloud.com (public.boxcloud.com)... 103.116.4.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|103.116.4.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 6778713048 (6.3G) [application/octet-stream]
    Saving to: ‘SRR6470907_S1_L001_R2_001.fastq.gz’
    
    SRR6470907_S1_L001_ 100%[===================>]   6.31G  14.5MB/s    in 7m 25s  
    
    2019-12-12 17:21:07 (14.5 MB/s) - ‘SRR6470907_S1_L001_R2_001.fastq.gz’ saved [6778713048/6778713048]
    
    --2019-12-12 17:21:09--  https://caltech.box.com/shared/static/ln14jjd4tz3hvgxf8zj2kmokof7f1nrf.gz
    Resolving caltech.box.com (caltech.box.com)... 103.116.4.197
    Connecting to caltech.box.com (caltech.box.com)|103.116.4.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/ln14jjd4tz3hvgxf8zj2kmokof7f1nrf.gz [following]
    --2019-12-12 17:21:10--  https://caltech.box.com/public/static/ln14jjd4tz3hvgxf8zj2kmokof7f1nrf.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/ln14jjd4tz3hvgxf8zj2kmokof7f1nrf.gz [following]
    --2019-12-12 17:21:10--  https://caltech.app.box.com/public/static/ln14jjd4tz3hvgxf8zj2kmokof7f1nrf.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 103.116.4.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|103.116.4.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!-CFZSuN_Qtpd3jFo3OMS6UdIb6n8hnhEU92hCLHKOzl5WL78fc6YEvnTsGBMC2gWBUaMdvkMq0UT3RkegIrKuohCJNeyRb6UXIzXAIPsuVVwmwkhcgPRIXn4KXm5Nub_VUGTAUtTEELuqbjlLurA1koY05AfijpRkdMrd6dq_e0L5SDjGIbn_d4LdtVY8Ty1Bs7m49XJDfc5xo74Myaomf-gm3EY02sfJ2JCJmcRMmZdDgoDaRUeAS8-w-lPZw7gLmWGw24gHpv2x0EdrSp_3mNEfPrTmZTHNUFDzCOtvRU4T0C3V1FU0lvBenkIs8wBLT9D4geMEjj6L_qYkbr6ZVVsq215oLGKsaZPWvvDi0ruJ2xE1FvO0k3cxG8K1xOesQkyb2Q4HMaueLISFYwjYtYoaptLI3og6YP-mSEBEmX8AB0vG2wtQTn3tU1hrMhvaHLa7z64RRbpyfr6r95EvViGiKCndPz20jZLcboZQ91nKX3A6lLkMaKKLG1aLQp6YcU1NUvb0616_X37Ij2THCsqK5VKVVOvZUuQfwSGyk6AKunp4L2T4awbfui96dFNXi_exzuC6CJdmAkOBLZeGMEXWGLWKJSZHCMcHjNArMJH24xb7JX_Muw1PJMbY7p8xzaJZ8us9MC5vEfC6bOvEub6EuegKZcbLNJsVirxbw7lGV3NdTxhZv-lTw0PgxCCf0Z4Obvy3ZesMmwp8GhSKT34t3FCDHMJhDdmN3LzaafCRt32mNvVIPi8-LZqdZI8U17nzAtqRvoax0b6nVX-7tLs8f6jlpfkWjCiUq6YaN5dl6cpFriTRRukWUkVCUbD7piKB9yFL0aIMvZcyjD_W-NAWUlsQxYWMeR-Yx_ZlBBRpmRR6Tu2xPQwQd9fbWUMeyB5C74o1V0lc4T_OGCU7MqE9Jo41btnvXbHgBSPPyGQ3zNHVGVYs_pCCVbedtjhaWmhHlztEcnRldao7h2z40gXeBuJAjK87MZ_ybBg_w_NIXr3E7iLG9HmkKiJngslzENTaFdR85PoqpPNroLt34YQQijizL1elyowlab3VnosScPMDKG5uvu6FYbiOQLIIK-znuzsMrej6M9RgNjv9HUW_aNv2-EQWhaOARHBt20p1CHBWseSeBdYKzTZurPZbS2HLTuzbGws4CuhcPvbaNJLyV9BTMNAG_rJqGT4gsisBJBWRnlpEJX0OahbDf2QPOobrNAj8rVQIoOeLXqqJzstgOyzx0gkfHszPx2fIQmRyRC-a2x0wrmCMUQWJ3kwHTKKVTSaM03_MVDpM7s5flxshlFqOpjWZsps235p6dqkd0-fKsL_4nNdym8sR4abwegEM1qIOTGDUO5cNw7UkKNFxuD3sZ39Si_5SEbZmpGgCgXfHlhYA9PV1g4qCLSDyVrpAh41Em_iTQouEb9JNcIwpmv11Tw./download [following]
    --2019-12-12 17:21:11--  https://public.boxcloud.com/d/1/b1!-CFZSuN_Qtpd3jFo3OMS6UdIb6n8hnhEU92hCLHKOzl5WL78fc6YEvnTsGBMC2gWBUaMdvkMq0UT3RkegIrKuohCJNeyRb6UXIzXAIPsuVVwmwkhcgPRIXn4KXm5Nub_VUGTAUtTEELuqbjlLurA1koY05AfijpRkdMrd6dq_e0L5SDjGIbn_d4LdtVY8Ty1Bs7m49XJDfc5xo74Myaomf-gm3EY02sfJ2JCJmcRMmZdDgoDaRUeAS8-w-lPZw7gLmWGw24gHpv2x0EdrSp_3mNEfPrTmZTHNUFDzCOtvRU4T0C3V1FU0lvBenkIs8wBLT9D4geMEjj6L_qYkbr6ZVVsq215oLGKsaZPWvvDi0ruJ2xE1FvO0k3cxG8K1xOesQkyb2Q4HMaueLISFYwjYtYoaptLI3og6YP-mSEBEmX8AB0vG2wtQTn3tU1hrMhvaHLa7z64RRbpyfr6r95EvViGiKCndPz20jZLcboZQ91nKX3A6lLkMaKKLG1aLQp6YcU1NUvb0616_X37Ij2THCsqK5VKVVOvZUuQfwSGyk6AKunp4L2T4awbfui96dFNXi_exzuC6CJdmAkOBLZeGMEXWGLWKJSZHCMcHjNArMJH24xb7JX_Muw1PJMbY7p8xzaJZ8us9MC5vEfC6bOvEub6EuegKZcbLNJsVirxbw7lGV3NdTxhZv-lTw0PgxCCf0Z4Obvy3ZesMmwp8GhSKT34t3FCDHMJhDdmN3LzaafCRt32mNvVIPi8-LZqdZI8U17nzAtqRvoax0b6nVX-7tLs8f6jlpfkWjCiUq6YaN5dl6cpFriTRRukWUkVCUbD7piKB9yFL0aIMvZcyjD_W-NAWUlsQxYWMeR-Yx_ZlBBRpmRR6Tu2xPQwQd9fbWUMeyB5C74o1V0lc4T_OGCU7MqE9Jo41btnvXbHgBSPPyGQ3zNHVGVYs_pCCVbedtjhaWmhHlztEcnRldao7h2z40gXeBuJAjK87MZ_ybBg_w_NIXr3E7iLG9HmkKiJngslzENTaFdR85PoqpPNroLt34YQQijizL1elyowlab3VnosScPMDKG5uvu6FYbiOQLIIK-znuzsMrej6M9RgNjv9HUW_aNv2-EQWhaOARHBt20p1CHBWseSeBdYKzTZurPZbS2HLTuzbGws4CuhcPvbaNJLyV9BTMNAG_rJqGT4gsisBJBWRnlpEJX0OahbDf2QPOobrNAj8rVQIoOeLXqqJzstgOyzx0gkfHszPx2fIQmRyRC-a2x0wrmCMUQWJ3kwHTKKVTSaM03_MVDpM7s5flxshlFqOpjWZsps235p6dqkd0-fKsL_4nNdym8sR4abwegEM1qIOTGDUO5cNw7UkKNFxuD3sZ39Si_5SEbZmpGgCgXfHlhYA9PV1g4qCLSDyVrpAh41Em_iTQouEb9JNcIwpmv11Tw./download
    Resolving public.boxcloud.com (public.boxcloud.com)... 103.116.4.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|103.116.4.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 4853988957 (4.5G) [application/octet-stream]
    Saving to: ‘SRR6470907_S1_L002_R1_001.fastq.gz’
    
    SRR6470907_S1_L002_ 100%[===================>]   4.52G  14.3MB/s    in 5m 28s  
    
    2019-12-12 17:26:39 (14.1 MB/s) - ‘SRR6470907_S1_L002_R1_001.fastq.gz’ saved [4853988957/4853988957]
    
    --2019-12-12 17:26:41--  https://caltech.box.com/shared/static/o5bwf9u2g7egi02by3e3hbvov8fgwbb3.gz
    Resolving caltech.box.com (caltech.box.com)... 103.116.4.197
    Connecting to caltech.box.com (caltech.box.com)|103.116.4.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/o5bwf9u2g7egi02by3e3hbvov8fgwbb3.gz [following]
    --2019-12-12 17:26:41--  https://caltech.box.com/public/static/o5bwf9u2g7egi02by3e3hbvov8fgwbb3.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/o5bwf9u2g7egi02by3e3hbvov8fgwbb3.gz [following]
    --2019-12-12 17:26:41--  https://caltech.app.box.com/public/static/o5bwf9u2g7egi02by3e3hbvov8fgwbb3.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 103.116.4.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|103.116.4.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!IVYOgKcGOzbilk9FlsGm88xV2io9D84Ow22xbB96OKsLZNYRxkZSs_jyMZQ8OC1FhBpCAB3dHUjCJCJ4NiaaUtOL8OwMjtQbMqL45-TTdPw8PpbCH4Y1hf9G1XrG336NJLhY7hwQEir1hiLmBGIuf4t9TrjHSkXx6GCbMdK9XYSafYznfKjR4usiOwhGcxPSTD96PW3BGPM8YPM9wfqeDJtWQymNJOGOGUiWaa6-gPrejYp1potoKKu-wWa6CMMxoO3aZ9pK2ZhLQzEpsbHmUkz2xB36RY5vqatTV1V58nxen06RjS1cExNfQ2HWxSS65IvfXmodDCiUEKeCUjHM08LmaPmvDdJfN6x9bIpGRyFWYqW845YcbvdWvvFA_TIfXcziixmBRiqc_kG4KqmIyV_A9hm40gANqeeb9UK-Oz41SX7kHRJp7Of_Th5SIctBrLN7NRdb-zqMmThHmnvCWJUk5eVSO0Xn5tOOJoiYLFl2PbiR3PYR0q-J6lkVUJvxU47-GBGWwd2S3yj115w7bXVfm77Ow895bSC4g6HdF21w2fw2CFhf5-XtDgp2W5TnxQjb_ZU0pNl7ZyJ8Vmf3b8eT1bPOfQT8UnMwZjGpcmXDsaRi7pxtzY7DP8nLVuY9HelV57_QIvEpOHK7hd810sfYu2r_kgVQRqrVaev9_LXZIh-aH2azE1YFjxUPrHwyETg4KLKVDpcQQC_Dj2NCWor2ThwUROC-O4_j5Zs4CYqMzPF6XiNBBOGu-zsjvDLTdOJ9ID8tqSGhsBo-Gkz5mWDz4st7A0kK4GjZ3jNnip02Q1Cf5s6AhDd8GdEi-45acxUTCB2rhcrxTTHV-XuCSTRQ7qQ0HNMMMFxk8CVNd6GNm_5_BaJqjd6NHpumIp61eBNfvdQIaU_ugiHNZiTnsVoCrTfk8Yn3ehSS0EZPHoSbvUwW9RFLP1sZk81836XQ3xHUiYjQ_owVntdqmBInhdq9lsxeB1acv5JLrpAPBpPSdOoAoKcSMFMEWMM36c5k9YNk7zj7Pi5h672oUzgP1mNj05hHyIgTQoXOjPtS9A6vcy8KDOOXSGQ7a9hXOeIw4nJ4U6ErlhZlvlRZ3RpJRsgavTXdAi6A5A-nrxoQVQvd3-1uwNQeEZJ5zGPsspbs51GZJOFt-kGVF2nqQ3K_w4vgMggST_V6Lgtr2wX0z1gMwyObrrCpgN0WJzSz4881mI9WNJedPc33Xiq7U33y02yRC1NowEQuIqJhPZZprH98e1GPgPnbaRQ-1im--XMpVaZQjWj71ptQvEvP_3b4JICWRIX9f4L_5PdEmfVkjrS4Khxd6HYekoiitQbm6d8auOdEp7PseYO-0yh0cwy1RCGVLKvrAp--bKhen4j4aVD5repULhU86R30xg5qJQWqp3Fvidw03bdao3tmtV7Pm9j9AIO3BaE./download [following]
    --2019-12-12 17:26:42--  https://public.boxcloud.com/d/1/b1!IVYOgKcGOzbilk9FlsGm88xV2io9D84Ow22xbB96OKsLZNYRxkZSs_jyMZQ8OC1FhBpCAB3dHUjCJCJ4NiaaUtOL8OwMjtQbMqL45-TTdPw8PpbCH4Y1hf9G1XrG336NJLhY7hwQEir1hiLmBGIuf4t9TrjHSkXx6GCbMdK9XYSafYznfKjR4usiOwhGcxPSTD96PW3BGPM8YPM9wfqeDJtWQymNJOGOGUiWaa6-gPrejYp1potoKKu-wWa6CMMxoO3aZ9pK2ZhLQzEpsbHmUkz2xB36RY5vqatTV1V58nxen06RjS1cExNfQ2HWxSS65IvfXmodDCiUEKeCUjHM08LmaPmvDdJfN6x9bIpGRyFWYqW845YcbvdWvvFA_TIfXcziixmBRiqc_kG4KqmIyV_A9hm40gANqeeb9UK-Oz41SX7kHRJp7Of_Th5SIctBrLN7NRdb-zqMmThHmnvCWJUk5eVSO0Xn5tOOJoiYLFl2PbiR3PYR0q-J6lkVUJvxU47-GBGWwd2S3yj115w7bXVfm77Ow895bSC4g6HdF21w2fw2CFhf5-XtDgp2W5TnxQjb_ZU0pNl7ZyJ8Vmf3b8eT1bPOfQT8UnMwZjGpcmXDsaRi7pxtzY7DP8nLVuY9HelV57_QIvEpOHK7hd810sfYu2r_kgVQRqrVaev9_LXZIh-aH2azE1YFjxUPrHwyETg4KLKVDpcQQC_Dj2NCWor2ThwUROC-O4_j5Zs4CYqMzPF6XiNBBOGu-zsjvDLTdOJ9ID8tqSGhsBo-Gkz5mWDz4st7A0kK4GjZ3jNnip02Q1Cf5s6AhDd8GdEi-45acxUTCB2rhcrxTTHV-XuCSTRQ7qQ0HNMMMFxk8CVNd6GNm_5_BaJqjd6NHpumIp61eBNfvdQIaU_ugiHNZiTnsVoCrTfk8Yn3ehSS0EZPHoSbvUwW9RFLP1sZk81836XQ3xHUiYjQ_owVntdqmBInhdq9lsxeB1acv5JLrpAPBpPSdOoAoKcSMFMEWMM36c5k9YNk7zj7Pi5h672oUzgP1mNj05hHyIgTQoXOjPtS9A6vcy8KDOOXSGQ7a9hXOeIw4nJ4U6ErlhZlvlRZ3RpJRsgavTXdAi6A5A-nrxoQVQvd3-1uwNQeEZJ5zGPsspbs51GZJOFt-kGVF2nqQ3K_w4vgMggST_V6Lgtr2wX0z1gMwyObrrCpgN0WJzSz4881mI9WNJedPc33Xiq7U33y02yRC1NowEQuIqJhPZZprH98e1GPgPnbaRQ-1im--XMpVaZQjWj71ptQvEvP_3b4JICWRIX9f4L_5PdEmfVkjrS4Khxd6HYekoiitQbm6d8auOdEp7PseYO-0yh0cwy1RCGVLKvrAp--bKhen4j4aVD5repULhU86R30xg5qJQWqp3Fvidw03bdao3tmtV7Pm9j9AIO3BaE./download
    Resolving public.boxcloud.com (public.boxcloud.com)... 103.116.4.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|103.116.4.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 8621336299 (8.0G) [application/octet-stream]
    Saving to: ‘SRR6470907_S1_L002_R2_001.fastq.gz’
    
    SRR6470907_S1_L002_ 100%[===================>]   8.03G  15.3MB/s    in 9m 1s   
    
    2019-12-12 17:35:44 (15.2 MB/s) - ‘SRR6470907_S1_L002_R2_001.fastq.gz’ saved [8621336299/8621336299]
    
    CPU times: user 17.6 s, sys: 3.11 s, total: 20.7 s
    Wall time: 50min 28s


Then, we verify the integrity of the files we downloaded to make sure they were not corrupted during the download.


```
!md5sum -c checksums.txt --ignore-missing
```

    SRR6470906_S1_L001_R1_001.fastq.gz: OK
    SRR6470906_S1_L001_R2_001.fastq.gz: OK
    SRR6470906_S1_L002_R1_001.fastq.gz: OK
    SRR6470906_S1_L002_R2_001.fastq.gz: OK
    SRR6470907_S1_L001_R1_001.fastq.gz: OK
    SRR6470907_S1_L001_R2_001.fastq.gz: OK
    SRR6470907_S1_L002_R1_001.fastq.gz: OK
    SRR6470907_S1_L002_R2_001.fastq.gz: OK


### Install `kb`

Install `kb` for running the kallisto|bustools workflow.


```
!pip install kb-python
```

    Collecting kb-python
    [?25l  Downloading https://files.pythonhosted.org/packages/62/c9/2e5b8fa2cd873a23ae1aeb128b33165d6a9387a2f56ea1fafec1d6d32477/kb_python-0.24.4-py3-none-any.whl (35.4MB)
    [K     |████████████████████████████████| 35.4MB 93.0MB/s 
    [?25hCollecting loompy>=3.0.6
    [?25l  Downloading https://files.pythonhosted.org/packages/36/52/74ed37ae5988522fbf87b856c67c4f80700e6452410b4cd80498c5f416f9/loompy-3.0.6.tar.gz (41kB)
    [K     |████████████████████████████████| 51kB 8.3MB/s 
    [?25hCollecting anndata>=0.6.22.post1
    [?25l  Downloading https://files.pythonhosted.org/packages/2b/72/87196c15f68d9865c31a43a10cf7c50bcbcedd5607d09f9aada0b3963103/anndata-0.6.22.post1-py3-none-any.whl (47kB)
    [K     |████████████████████████████████| 51kB 8.2MB/s 
    [?25hRequirement already satisfied: h5py in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (2.8.0)
    Requirement already satisfied: numpy in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (1.17.4)
    Requirement already satisfied: scipy in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (1.3.3)
    Requirement already satisfied: setuptools in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (42.0.2)
    Requirement already satisfied: numba in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (0.40.1)
    Requirement already satisfied: click in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (7.0)
    Collecting numpy-groupies
    [?25l  Downloading https://files.pythonhosted.org/packages/96/7a/2196465530e72084c6bb97cd49bf8ccdc83919cc94755727aa148effbc0f/numpy_groupies-0.9.9.tar.gz (43kB)
    [K     |████████████████████████████████| 51kB 7.3MB/s 
    [?25hRequirement already satisfied: pandas>=0.23.0 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python) (0.25.3)
    Requirement already satisfied: natsort in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python) (5.5.0)
    Requirement already satisfied: six in /usr/local/lib/python3.6/dist-packages (from h5py->loompy>=3.0.6->kb-python) (1.12.0)
    Requirement already satisfied: llvmlite>=0.25.0dev0 in /usr/local/lib/python3.6/dist-packages (from numba->loompy>=3.0.6->kb-python) (0.30.0)
    Requirement already satisfied: python-dateutil>=2.6.1 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python) (2.6.1)
    Requirement already satisfied: pytz>=2017.2 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python) (2018.9)
    Building wheels for collected packages: loompy, numpy-groupies
      Building wheel for loompy (setup.py) ... [?25l[?25hdone
      Created wheel for loompy: filename=loompy-3.0.6-cp36-none-any.whl size=47896 sha256=59c25a4deb5ede0f81a8e517a49472995f8e78dfba775923ecd967c87fe50386
      Stored in directory: /root/.cache/pip/wheels/f9/a4/90/5a98ad83419732b0fba533b81a2a52ba3dbe230a936ca4cdc9
      Building wheel for numpy-groupies (setup.py) ... [?25l[?25hdone
      Created wheel for numpy-groupies: filename=numpy_groupies-0+unknown-cp36-none-any.whl size=27879 sha256=87270c8976a5752b0f7ea76542caffdf0070a3f736e8220c086a16e2343863ec
      Stored in directory: /root/.cache/pip/wheels/6c/fb/3d/5c43eb691bd92a3ddd0ebeb6e7e78ceaf3ae1cb8d54b89a7fb
    Successfully built loompy numpy-groupies
    Installing collected packages: numpy-groupies, loompy, anndata, kb-python
    Successfully installed anndata-0.6.22.post1 kb-python-0.24.4 loompy-3.0.6 numpy-groupies-0+unknown


### Download a pre-built human RNA velocity index

`kb` provides a pre-built human RNA velocity index, which is suitable for 10x data only. See this notebook to understand how the index was generated: https://github.com/linnarsson-lab/loompy/blob/master/notebooks/build_index.ipynb Because this is a velocity index, we need to provide the `-c1` and `-c2` options when downloading to indicate what we want to name our transcripts-to-capture lists. These lists contain transcript IDs that correspond to spliced and unspliced variants of genes, the sequences of which were used to generate the `index.idx` file.

__Note:__ See [this notebook](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_velocity_index.ipynb) for a tutorial on how to build custom transcriptome or RNA velocity indices.


```
%%time
!kb ref -d linnarsson -i index.idx -g t2g.txt -c1 spliced_t2c.txt -c2 unspliced_t2c.txt
```

    [2019-12-12 17:40:34,599]    INFO Downloading files for linnarsson from https://caltech.box.com/shared/static/kyf7ai5s8y2l0vycl5yxunrappvrf0yx.gz to tmp/kyf7ai5s8y2l0vycl5yxunrappvrf0yx.gz
    [2019-12-12 17:45:54,084]    INFO Extracting files from tmp/kyf7ai5s8y2l0vycl5yxunrappvrf0yx.gz
    CPU times: user 1.32 s, sys: 174 ms, total: 1.5 s
    Wall time: 6min 30s


### Generate RNA velocity count matrices

The following command will generate an RNA count matrix of cells (rows) by genes (columns) in H5AD format, which is a binary format used to store [Anndata](https://anndata.readthedocs.io/en/stable/) objects. Notice we are providing the index and transcript-to-gene mapping we downloaded in the previous step to the `-i` and `-g` arguments respectively, as well as the transcripts-to-capture lists to the `-c1` and `-c2` arguments. Also, these reads were generated with the 10x Genomics Chromium Single Cell v2 Chemistry, hence the `-x 10xv2` argument. To view other supported technologies, run `kb --list`.

The `--filter` flag is used to filter out barcodes with low UMI counts. This will generate two matrices, one in the `counts_unfiltered` directory and another in the `counts_filtered` directory.

__Note:__ If you would like a Loom file instead, replace the `--h5ad` flag with `--loom`. If you want to use the raw matrix output by `kb` instead of their H5AD or Loom converted files, omit these flags.

#### SRR6470906


```
%%time
# Uncomment the following two lines if you want to use checkpoint files
# !wget https://caltech.box.com/shared/static/wt1xo5ro9v7qo4fi4dl52kbbcg0omlve.gz -O SRR64070906.tar.gz
# !tar -xzvf SRR64070906.tar.gz
!kb count --h5ad -i index.idx -g t2g.txt -x 10xv2 -o SRR6470906 \
-c1 spliced_t2c.txt -c2 unspliced_t2c.txt --lamanno --filter bustools -t 2 \
SRR6470906_S1_L001_R1_001.fastq.gz \
SRR6470906_S1_L001_R2_001.fastq.gz \
SRR6470906_S1_L002_R1_001.fastq.gz \
SRR6470906_S1_L002_R2_001.fastq.gz
```

    [2019-12-12 17:59:33,838]    INFO Generating BUS file from
    [2019-12-12 17:59:33,838]    INFO         SRR6470906_S1_L001_R1_001.fastq.gz
    [2019-12-12 17:59:33,838]    INFO         SRR6470906_S1_L001_R2_001.fastq.gz
    [2019-12-12 17:59:33,838]    INFO         SRR6470906_S1_L002_R1_001.fastq.gz
    [2019-12-12 17:59:33,838]    INFO         SRR6470906_S1_L002_R2_001.fastq.gz
    [2019-12-12 18:29:18,309]    INFO Sorting BUS file SRR6470906/output.bus to tmp/output.s.bus
    [2019-12-12 18:30:52,034]    INFO Whitelist not provided
    [2019-12-12 18:30:52,034]    INFO Copying pre-packaged 10XV2 whitelist to SRR6470906
    [2019-12-12 18:30:52,217]    INFO Inspecting BUS file tmp/output.s.bus
    [2019-12-12 18:31:19,945]    INFO Correcting BUS records in tmp/output.s.bus to tmp/output.s.c.bus with whitelist SRR6470906/10xv2_whitelist.txt
    [2019-12-12 18:31:44,539]    INFO Sorting BUS file tmp/output.s.c.bus to SRR6470906/output.unfiltered.bus
    [2019-12-12 18:32:15,671]    INFO Capturing records from BUS file SRR6470906/output.unfiltered.bus to tmp/spliced.bus with capture list spliced_t2c.txt
    [2019-12-12 18:32:55,811]    INFO Sorting BUS file tmp/spliced.bus to SRR6470906/spliced.unfiltered.bus
    [2019-12-12 18:33:18,574]    INFO Generating count matrix SRR6470906/counts_unfiltered/spliced from BUS file SRR6470906/spliced.unfiltered.bus
    [2019-12-12 18:33:47,704]    INFO Capturing records from BUS file SRR6470906/output.unfiltered.bus to tmp/unspliced.bus with capture list unspliced_t2c.txt
    [2019-12-12 18:34:19,695]    INFO Sorting BUS file tmp/unspliced.bus to SRR6470906/unspliced.unfiltered.bus
    [2019-12-12 18:34:30,651]    INFO Generating count matrix SRR6470906/counts_unfiltered/unspliced from BUS file SRR6470906/unspliced.unfiltered.bus
    [2019-12-12 18:35:28,217]    INFO Writing matrices to h5ad SRR6470906/counts_unfiltered/adata.h5ad
    [2019-12-12 18:35:28,762]    INFO Filtering with bustools
    [2019-12-12 18:35:28,763]    INFO Generating whitelist SRR6470906/filter_barcodes.txt from BUS file SRR6470906/output.unfiltered.bus
    [2019-12-12 18:35:29,383]    INFO Capturing records from BUS file SRR6470906/output.unfiltered.bus to tmp/output.filtered.bus with capture list SRR6470906/filter_barcodes.txt
    [2019-12-12 18:35:37,452]    INFO Sorting BUS file tmp/output.filtered.bus to SRR6470906/output.filtered.bus
    [2019-12-12 18:36:01,835]    INFO Capturing records from BUS file SRR6470906/output.filtered.bus to tmp/spliced.bus with capture list spliced_t2c.txt
    [2019-12-12 18:36:38,648]    INFO Sorting BUS file tmp/spliced.bus to SRR6470906/spliced.filtered.bus
    [2019-12-12 18:36:57,294]    INFO Generating count matrix SRR6470906/counts_filtered/spliced from BUS file SRR6470906/spliced.filtered.bus
    [2019-12-12 18:37:23,234]    INFO Capturing records from BUS file SRR6470906/output.filtered.bus to tmp/unspliced.bus with capture list unspliced_t2c.txt
    [2019-12-12 18:37:54,444]    INFO Sorting BUS file tmp/unspliced.bus to SRR6470906/unspliced.filtered.bus
    [2019-12-12 18:38:04,159]    INFO Generating count matrix SRR6470906/counts_filtered/unspliced from BUS file SRR6470906/unspliced.filtered.bus
    [2019-12-12 18:38:51,739]    INFO Writing matrices to h5ad SRR6470906/counts_filtered/adata.h5ad
    CPU times: user 8.86 s, sys: 1.1 s, total: 9.97 s
    Wall time: 39min 20s


#### SRR6470907


```
%%time
# Uncomment the following two lines if you want to use checkpoint files
# !wget https://caltech.box.com/shared/static/bulmquvaflr18v0nqb8c3ataf9di55r8.gz -O SRR64070907.tar.gz
# !tar -xzvf SRR64070907.tar.gz
!kb count --h5ad -i index.idx -g t2g.txt -x 10xv2 -o SRR6470907 \
-c1 spliced_t2c.txt -c2 unspliced_t2c.txt --lamanno --filter bustools \
SRR6470907_S1_L001_R1_001.fastq.gz \
SRR6470907_S1_L001_R2_001.fastq.gz \
SRR6470907_S1_L002_R1_001.fastq.gz \
SRR6470907_S1_L002_R2_001.fastq.gz
```

    [2019-12-12 18:38:55,843]    INFO Generating BUS file from
    [2019-12-12 18:38:55,843]    INFO         SRR6470907_S1_L001_R1_001.fastq.gz
    [2019-12-12 18:38:55,844]    INFO         SRR6470907_S1_L001_R2_001.fastq.gz
    [2019-12-12 18:38:55,844]    INFO         SRR6470907_S1_L002_R1_001.fastq.gz
    [2019-12-12 18:38:55,844]    INFO         SRR6470907_S1_L002_R2_001.fastq.gz
    [2019-12-12 19:10:35,733]    INFO Sorting BUS file SRR6470907/output.bus to tmp/output.s.bus
    [2019-12-12 19:12:20,479]    INFO Whitelist not provided
    [2019-12-12 19:12:20,479]    INFO Copying pre-packaged 10XV2 whitelist to SRR6470907
    [2019-12-12 19:12:20,592]    INFO Inspecting BUS file tmp/output.s.bus
    [2019-12-12 19:12:54,911]    INFO Correcting BUS records in tmp/output.s.bus to tmp/output.s.c.bus with whitelist SRR6470907/10xv2_whitelist.txt
    [2019-12-12 19:13:22,556]    INFO Sorting BUS file tmp/output.s.c.bus to SRR6470907/output.unfiltered.bus
    [2019-12-12 19:14:04,585]    INFO Capturing records from BUS file SRR6470907/output.unfiltered.bus to tmp/spliced.bus with capture list spliced_t2c.txt
    [2019-12-12 19:14:52,037]    INFO Sorting BUS file tmp/spliced.bus to SRR6470907/spliced.unfiltered.bus
    [2019-12-12 19:15:19,416]    INFO Generating count matrix SRR6470907/counts_unfiltered/spliced from BUS file SRR6470907/spliced.unfiltered.bus
    [2019-12-12 19:15:51,857]    INFO Capturing records from BUS file SRR6470907/output.unfiltered.bus to tmp/unspliced.bus with capture list unspliced_t2c.txt
    [2019-12-12 19:16:30,927]    INFO Sorting BUS file tmp/unspliced.bus to SRR6470907/unspliced.unfiltered.bus
    [2019-12-12 19:16:46,008]    INFO Generating count matrix SRR6470907/counts_unfiltered/unspliced from BUS file SRR6470907/unspliced.unfiltered.bus
    [2019-12-12 19:17:57,444]    INFO Writing matrices to h5ad SRR6470907/counts_unfiltered/adata.h5ad
    [2019-12-12 19:17:58,086]    INFO Filtering with bustools
    [2019-12-12 19:17:58,086]    INFO Generating whitelist SRR6470907/filter_barcodes.txt from BUS file SRR6470907/output.unfiltered.bus
    [2019-12-12 19:17:58,897]    INFO Capturing records from BUS file SRR6470907/output.unfiltered.bus to tmp/output.filtered.bus with capture list SRR6470907/filter_barcodes.txt
    [2019-12-12 19:18:08,972]    INFO Sorting BUS file tmp/output.filtered.bus to SRR6470907/output.filtered.bus
    [2019-12-12 19:18:42,842]    INFO Capturing records from BUS file SRR6470907/output.filtered.bus to tmp/spliced.bus with capture list spliced_t2c.txt
    [2019-12-12 19:19:26,502]    INFO Sorting BUS file tmp/spliced.bus to SRR6470907/spliced.filtered.bus
    [2019-12-12 19:19:50,027]    INFO Generating count matrix SRR6470907/counts_filtered/spliced from BUS file SRR6470907/spliced.filtered.bus
    [2019-12-12 19:20:19,058]    INFO Capturing records from BUS file SRR6470907/output.filtered.bus to tmp/unspliced.bus with capture list unspliced_t2c.txt
    [2019-12-12 19:20:56,254]    INFO Sorting BUS file tmp/unspliced.bus to SRR6470907/unspliced.filtered.bus
    [2019-12-12 19:21:08,867]    INFO Generating count matrix SRR6470907/counts_filtered/unspliced from BUS file SRR6470907/unspliced.filtered.bus
    [2019-12-12 19:22:11,325]    INFO Writing matrices to h5ad SRR6470907/counts_filtered/adata.h5ad
    CPU times: user 9.69 s, sys: 1.25 s, total: 10.9 s
    Wall time: 43min 24s


## Analysis

In this part of the tutorial, we will load the RNA count matrix generated by `kb count` into Python and cluster the cells with Louvain.

### Install packages

Google Colab does not come with `scanpy`, `python-igraph`, or `louvain` (but comes with `matplotlib`, `numpy`, `pandas`, and `scipy`).


```
!pip install scanpy loompy scvelo anndata velocyto python-igraph louvain
```

    Collecting scanpy
    [?25l  Downloading https://files.pythonhosted.org/packages/c9/df/e0957422a85e924501afbff2f9283e02d93cefcb04aa8319ade7af68d266/scanpy-1.4.4.post1-py3-none-any.whl (1.9MB)
    [K     |████████████████████████████████| 1.9MB 2.7MB/s 
    [?25hRequirement already satisfied: loompy in /usr/local/lib/python3.6/dist-packages (3.0.6)
    Collecting scvelo
    [?25l  Downloading https://files.pythonhosted.org/packages/62/da/0b40b56a7a7c2b01fb3d978b7225c3f33558092a8c7a8374d67285f04c4b/scvelo-0.1.24-py3-none-any.whl (168kB)
    [K     |████████████████████████████████| 174kB 71.2MB/s 
    [?25hRequirement already satisfied: anndata in /usr/local/lib/python3.6/dist-packages (0.6.22.post1)
    Collecting velocyto
    [?25l  Downloading https://files.pythonhosted.org/packages/81/66/e8fff9d3b824fd99c0f678c47c740fec058ce2d9a0cfdf11b114ea8889f2/velocyto-0.17.17.tar.gz (198kB)
    [K     |████████████████████████████████| 204kB 70.5MB/s 
    [?25hRequirement already satisfied: scikit-learn!=0.21.0,!=0.21.1,>=0.19.1 in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.21.3)
    Requirement already satisfied: scipy>=1.3 in /usr/local/lib/python3.6/dist-packages (from scanpy) (1.3.3)
    Requirement already satisfied: natsort in /usr/local/lib/python3.6/dist-packages (from scanpy) (5.5.0)
    Requirement already satisfied: pandas>=0.21 in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.25.3)
    Requirement already satisfied: statsmodels>=0.10.0rc2 in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.10.2)
    Requirement already satisfied: seaborn in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.9.0)
    Requirement already satisfied: patsy in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.5.1)
    Collecting matplotlib==3.0.*
    [?25l  Downloading https://files.pythonhosted.org/packages/e9/69/f5e05f578585ed9935247be3788b374f90701296a70c8871bcd6d21edb00/matplotlib-3.0.3-cp36-cp36m-manylinux1_x86_64.whl (13.0MB)
    [K     |████████████████████████████████| 13.0MB 47.2MB/s 
    [?25hRequirement already satisfied: networkx in /usr/local/lib/python3.6/dist-packages (from scanpy) (2.4)
    Requirement already satisfied: tqdm in /usr/local/lib/python3.6/dist-packages (from scanpy) (4.28.1)
    Requirement already satisfied: joblib in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.14.0)
    Collecting numba>=0.41.0
    [?25l  Downloading https://files.pythonhosted.org/packages/53/34/22b6c2ded558976b5367be01b851ae679a0d1ba994de002d54afe84187b5/numba-0.46.0-cp36-cp36m-manylinux1_x86_64.whl (3.6MB)
    [K     |████████████████████████████████| 3.6MB 47.5MB/s 
    [?25hRequirement already satisfied: umap-learn>=0.3.0 in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.3.10)
    Requirement already satisfied: importlib-metadata>=0.7; python_version < "3.8" in /usr/local/lib/python3.6/dist-packages (from scanpy) (1.2.0)
    Requirement already satisfied: tables in /usr/local/lib/python3.6/dist-packages (from scanpy) (3.4.4)
    Requirement already satisfied: h5py in /usr/local/lib/python3.6/dist-packages (from scanpy) (2.8.0)
    Requirement already satisfied: setuptools in /usr/local/lib/python3.6/dist-packages (from loompy) (42.0.2)
    Requirement already satisfied: click in /usr/local/lib/python3.6/dist-packages (from loompy) (7.0)
    Requirement already satisfied: numpy-groupies in /usr/local/lib/python3.6/dist-packages (from loompy) (0+unknown)
    Requirement already satisfied: numpy in /usr/local/lib/python3.6/dist-packages (from loompy) (1.17.4)
    Requirement already satisfied: cython in /usr/local/lib/python3.6/dist-packages (from velocyto) (0.29.14)
    Collecting pysam
    [?25l  Downloading https://files.pythonhosted.org/packages/15/e7/2dab8bb0ac739555e69586f1492f0ff6bc4a1f8312992a83001d3deb77ac/pysam-0.15.3.tar.gz (3.2MB)
    [K     |████████████████████████████████| 3.2MB 39.3MB/s 
    [?25hRequirement already satisfied: pytz>=2017.2 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.21->scanpy) (2018.9)
    Requirement already satisfied: python-dateutil>=2.6.1 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.21->scanpy) (2.6.1)
    Requirement already satisfied: six in /usr/local/lib/python3.6/dist-packages (from patsy->scanpy) (1.12.0)
    Requirement already satisfied: pyparsing!=2.0.4,!=2.1.2,!=2.1.6,>=2.0.1 in /usr/local/lib/python3.6/dist-packages (from matplotlib==3.0.*->scanpy) (2.4.5)
    Requirement already satisfied: cycler>=0.10 in /usr/local/lib/python3.6/dist-packages (from matplotlib==3.0.*->scanpy) (0.10.0)
    Requirement already satisfied: kiwisolver>=1.0.1 in /usr/local/lib/python3.6/dist-packages (from matplotlib==3.0.*->scanpy) (1.1.0)
    Requirement already satisfied: decorator>=4.3.0 in /usr/local/lib/python3.6/dist-packages (from networkx->scanpy) (4.4.1)
    Requirement already satisfied: llvmlite>=0.30.0dev0 in /usr/local/lib/python3.6/dist-packages (from numba>=0.41.0->scanpy) (0.30.0)
    Requirement already satisfied: zipp>=0.5 in /usr/local/lib/python3.6/dist-packages (from importlib-metadata>=0.7; python_version < "3.8"->scanpy) (0.6.0)
    Requirement already satisfied: numexpr>=2.5.2 in /usr/local/lib/python3.6/dist-packages (from tables->scanpy) (2.7.0)
    Requirement already satisfied: more-itertools in /usr/local/lib/python3.6/dist-packages (from zipp>=0.5->importlib-metadata>=0.7; python_version < "3.8"->scanpy) (8.0.0)
    Building wheels for collected packages: velocyto, pysam
      Building wheel for velocyto (setup.py) ... [?25l[?25hdone
      Created wheel for velocyto: filename=velocyto-0.17.17-cp36-cp36m-linux_x86_64.whl size=367926 sha256=32f44ce715f59ffb9f36c324fc2cd0472dcd7507754fb1f8ef8868bf8c8b38c6
      Stored in directory: /root/.cache/pip/wheels/e0/10/47/5a2aa6a7179b17b50a19cdba1df71798ade77e7d9ce98c5300
      Building wheel for pysam (setup.py) ... [?25l[?25hdone
      Created wheel for pysam: filename=pysam-0.15.3-cp36-cp36m-linux_x86_64.whl size=8790482 sha256=3abfe873fccbf84384bed55c40e200be7425d00cfaa6c5877eddc1ca248add6b
      Stored in directory: /root/.cache/pip/wheels/85/ab/84/86ca6dda37a6fc85687b67be7345b735cd82f6584bea56f327
    Successfully built velocyto pysam
    [31mERROR: albumentations 0.1.12 has requirement imgaug<0.2.7,>=0.2.5, but you'll have imgaug 0.2.9 which is incompatible.[0m
    Installing collected packages: matplotlib, numba, scanpy, scvelo, pysam, velocyto
      Found existing installation: matplotlib 3.1.2
        Uninstalling matplotlib-3.1.2:
          Successfully uninstalled matplotlib-3.1.2
      Found existing installation: numba 0.40.1
        Uninstalling numba-0.40.1:
          Successfully uninstalled numba-0.40.1
    Successfully installed matplotlib-3.0.3 numba-0.46.0 pysam-0.15.3 scanpy-1.4.4.post1 scvelo-0.1.24 velocyto-0.17.17




We also install an `R` package called `princurve`, which is used to fit a principal curve.

__Note:__ Google Colab does not official support `R` as of writing this tutorial, so this process may change with the official release. At the moment, we need to "hack" Google Colab to install `R` dependencies.


```
!Rscript -e "install.packages('princurve')"
```

    Installing package into ‘/usr/local/lib/R/site-library’
    (as ‘lib’ is unspecified)
    trying URL 'https://cran.rstudio.com/src/contrib/princurve_2.1.4.tar.gz'
    Content type 'application/x-gzip' length 42249 bytes (41 KB)
    ==================================================
    downloaded 41 KB
    
    * installing *source* package ‘princurve’ ...
    ** package ‘princurve’ successfully unpacked and MD5 sums checked
    ** using staged installation
    ** libs
    g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG  -I"/usr/lib/R/site-library/Rcpp/include"   -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-uuRxut/r-base-3.6.1=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c RcppExports.cpp -o RcppExports.o
    g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG  -I"/usr/lib/R/site-library/Rcpp/include"   -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-uuRxut/r-base-3.6.1=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c project_to_curve.cpp -o project_to_curve.o
    g++ -std=gnu++11 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions -Wl,-z,relro -o princurve.so RcppExports.o project_to_curve.o -L/usr/lib/R/lib -lR
    installing to /usr/local/lib/R/site-library/00LOCK-princurve/00new/princurve/libs
    ** R
    ** inst
    ** byte-compile and prepare package for lazy loading
    ** help
    *** installing help indices
    *** copying figures
    ** building package indices
    ** testing if installed package can be loaded from temporary location
    ** checking absolute paths in shared objects and dynamic libraries
    ** testing if installed package can be loaded from final location
    ** testing if installed package keeps a record of temporary installation path
    * DONE (princurve)
    
    The downloaded source packages are in
    	‘/tmp/RtmpNrrQA5/downloaded_packages’


### Import packages


```
import anndata
import igraph
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import rpy2.robjects as robj
import scanpy as sc
import scipy
import scipy.optimize
import scvelo as scv
import sklearn
import velocyto as vcy

from collections import Counter
from IPython.core.display import display, HTML
from numpy_groupies import aggregate, aggregate_np
from rpy2.robjects.packages import importr
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import pdist, squareform

matplotlib.rcParams.update({'font.size': 12})
%config InlineBackend.figure_format = 'retina'
```

### Import H5AD-formatted Anndata matrices


```
adata_06 = anndata.read('SRR6470906/counts_unfiltered/adata.h5ad')
adata_07 = anndata.read('SRR6470907/counts_unfiltered/adata.h5ad')
```

### Combine both matrices into a single Anndata object


```
# Before we do, we need to make sure we can trace each row back to its
# original anndata object.
adata_06.obs['run'] = '06'
adata_07.obs['run'] = '07'
adata_06.obs['bcs'] = adata_06.obs.index
adata_07.obs['bcs'] = adata_07.obs.index
adata_06.obs.index = adata_06.obs['bcs'] + '.' + adata_06.obs['run']
adata_07.obs.index = adata_07.obs['bcs'] + '.' + adata_07.obs['run']

adata = adata_06.concatenate(adata_07, batch_key='batch')
```

Let's remove the `-0` prefixes from the `obs` index.


```
adata.obs.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>batch</th>
      <th>bcs</th>
      <th>run</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>AAACCTGAGAAACCGC.06-0</th>
      <td>0</td>
      <td>AAACCTGAGAAACCGC</td>
      <td>06</td>
    </tr>
    <tr>
      <th>AAACCTGAGAAGGACA.06-0</th>
      <td>0</td>
      <td>AAACCTGAGAAGGACA</td>
      <td>06</td>
    </tr>
    <tr>
      <th>AAACCTGAGAAGGCCT.06-0</th>
      <td>0</td>
      <td>AAACCTGAGAAGGCCT</td>
      <td>06</td>
    </tr>
    <tr>
      <th>AAACCTGAGAAGGGTA.06-0</th>
      <td>0</td>
      <td>AAACCTGAGAAGGGTA</td>
      <td>06</td>
    </tr>
    <tr>
      <th>AAACCTGAGAAGGTTT.06-0</th>
      <td>0</td>
      <td>AAACCTGAGAAGGTTT</td>
      <td>06</td>
    </tr>
  </tbody>
</table>
</div>




```
adata.obs.index = adata.obs.index.str.split('-').str[0]
adata.obs.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>batch</th>
      <th>bcs</th>
      <th>run</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>AAACCTGAGAAACCGC.06</th>
      <td>0</td>
      <td>AAACCTGAGAAACCGC</td>
      <td>06</td>
    </tr>
    <tr>
      <th>AAACCTGAGAAGGACA.06</th>
      <td>0</td>
      <td>AAACCTGAGAAGGACA</td>
      <td>06</td>
    </tr>
    <tr>
      <th>AAACCTGAGAAGGCCT.06</th>
      <td>0</td>
      <td>AAACCTGAGAAGGCCT</td>
      <td>06</td>
    </tr>
    <tr>
      <th>AAACCTGAGAAGGGTA.06</th>
      <td>0</td>
      <td>AAACCTGAGAAGGGTA</td>
      <td>06</td>
    </tr>
    <tr>
      <th>AAACCTGAGAAGGTTT.06</th>
      <td>0</td>
      <td>AAACCTGAGAAGGTTT</td>
      <td>06</td>
    </tr>
  </tbody>
</table>
</div>



Let's also remove the gene versions from the `var` index.


```
adata.var.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENSG00000277400.1</th>
    </tr>
    <tr>
      <th>ENSG00000274847.1</th>
    </tr>
    <tr>
      <th>ENSG00000276256.1</th>
    </tr>
    <tr>
      <th>ENSG00000278198.1</th>
    </tr>
    <tr>
      <th>ENSG00000273496.1</th>
    </tr>
  </tbody>
</table>
</div>




```
adata.var.index = adata.var.index.str.split('.').str[0]
adata.var.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENSG00000277400</th>
    </tr>
    <tr>
      <th>ENSG00000274847</th>
    </tr>
    <tr>
      <th>ENSG00000276256</th>
    </tr>
    <tr>
      <th>ENSG00000278198</th>
    </tr>
    <tr>
      <th>ENSG00000273496</th>
    </tr>
  </tbody>
</table>
</div>




```
adata
```




    AnnData object with n_obs × n_vars = 374814 × 58367 
        obs: 'batch', 'bcs', 'run'
        layers: 'spliced', 'unspliced'



### Select the right barcodes and genes

Our Anndata matrix so far contains all genes and barcodes. Here, we will take the genes and barcodes present in our matrix and the matrix published with this dataset.

#### Download the matrix


```
them = scv.read('data/hgForebrainGlut.loom', cleanup=True, sparse=True, cache=True, backup_url='http://pklab.med.harvard.edu/velocyto/hgForebrainGlut/hgForebrainGlut.loom')
them.var_names_make_unique()
```


    HBox(children=(IntProgress(value=1, bar_style='info', description='hgForebrainGlut.loom', max=1, style=Progres…


    


    Variable names are not unique. To make them unique, call `.var_names_make_unique`.
    ... storing 'Chromosome' as categorical
    ... storing 'Strand' as categorical


Clean it up so that we can compare it to ours


```
them.obs["bcs"] = them.obs.index.str.slice(11,-1)
them.obs["bid"] = them.obs.index.str.slice(8,10)
them.obs["run"] = them.obs.bid.map(lambda x: "06" if x=="29" else "07")
them.obs.index = them.obs.bcs.values + "."+ them.obs["run"]
```

#### Take the intersection


```
final = adata[adata.obs.index.isin(them.obs.index),:]
final = final[:,final.var.index.isin(them.var.Accession)]
```


```
final
```




    View of AnnData object with n_obs × n_vars = 1720 × 31029 
        obs: 'batch', 'bcs', 'run'
        layers: 'spliced', 'unspliced'



We also add an additional `ambiguous` layer, `CellID` and `Clusters` obs and `Accession` var, all of which are required by `Velocyto`.


```
final.layers["ambiguous"] = scipy.sparse.csr_matrix(np.zeros(final.X.shape))
final.obs["CellID"] = final.obs.index
final.obs["Clusters"] = final.obs.index.map(them.obs["Clusters"])
final.var["Accession"] = final.var.index
```


```
scv.pp.show_proportions(final)
print(final)
```

    Abundance of ['spliced', 'unspliced']: [0.61 0.39]
    AnnData object with n_obs × n_vars = 1720 × 31029 
        obs: 'batch', 'bcs', 'run', 'CellID', 'Clusters'
        var: 'Accession'
        layers: 'spliced', 'unspliced'


### Save the final Anndata as a Loom file


```
final.write_loom('final.loom')
```

    CPU times: user 4.74 s, sys: 237 ms, total: 4.98 s
    Wall time: 4.98 s


### Run the `Velocyto` python pipeline


```
# Some helper functions
def array_to_rmatrix(X):
    nr, nc = X.shape
    xvec = robj.FloatVector(X.transpose().reshape((X.size)))
    xr = robj.r.matrix(xvec, nrow=nr, ncol=nc)
    return xr

def principal_curve(X, pca=True):
    """
    input : numpy.array
    returns:
    Result::Object
        Methods:
        projections - the matrix of the projectiond
        ixsort - the order ot the points (as in argsort)
        arclength - the lenght of the arc from the beginning to the point
    """
    # convert array to R matrix
    xr = array_to_rmatrix(X)
    
    if pca:
        #perform pca
        t = robj.r.prcomp(xr)
        #determine dimensionality reduction
        usedcomp = max( sum( np.array(t[t.names.index('sdev')]) > 1.1) , 4)
        usedcomp = min([usedcomp, sum( np.array(t[t.names.index('sdev')]) > 0.25), X.shape[0]])
        Xpc = np.array(t[t.names.index('x')])[:,:usedcomp]
        # convert array to R matrix
        xr = array_to_rmatrix(Xpc)

    #import the correct namespace
    princurve = importr("princurve", on_conflict="warn")
    
    #call the function
    fit1 = princurve.principal_curve(xr)
    
    #extract the outputs
    class Results:
        pass
    results = Results()
    results.projections = np.array( fit1[0] )
    results.ixsort = np.array( fit1[1] ) - 1 # R is 1 indexed
    results.arclength = np.array( fit1[2] )
    results.dist = np.array( fit1[3] )
    
    if pca:
        results.PCs = np.array(xr) #only the used components
        
    return results
```

### Read the saved Loom file with `Velocyto`


```
vlm = vcy.VelocytoLoom("final.loom")
```

The rest of the notebook comes from the following notebook: https://github.com/velocyto-team/velocyto-notebooks/blob/master/python/hgForebrainGlutamatergic.ipynb


```
# Load an initial clustering (Louvein)
labels = vlm.ca["Clusters"]
manual_annotation = {str(i):[i] for i in labels}
annotation_dict = {v:k for k, values in manual_annotation.items() for v in values }
clusters = np.array([annotation_dict[i] for i in labels])
colors20 = np.vstack((plt.cm.tab20b(np.linspace(0., 1, 20))[::2], plt.cm.tab20c(np.linspace(0, 1, 20))[1::2]))
vlm.set_clusters(clusters, cluster_colors_dict={k:colors20[v[0] % 20,:] for k,v in manual_annotation.items()})
```


```
# just to find the initial cell size
vlm.normalize("S", size=True, log=False)
vlm.normalize("U", size=True,  log=False)
```


```
vlm.score_detection_levels(min_expr_counts=30, min_cells_express=20,
                           min_expr_counts_U=0, min_cells_express_U=0)
vlm.filter_genes(by_detection_levels=True)
```


```
vlm.score_cv_vs_mean(2000, plot=True, max_expr_avg=50, winsorize=True, winsor_perc=(1,99.8), svr_gamma=0.01, min_expr_cells=50)
vlm.filter_genes(by_cv_vs_mean=True)
```


![png](kb_velocity_files/kb_velocity_55_0.png)



```
vlm.score_detection_levels(min_expr_counts=0, min_cells_express=0,
                           min_expr_counts_U=25, min_cells_express_U=20)
vlm.score_cluster_expression(min_avg_U=0.007, min_avg_S=0.06)
vlm.filter_genes(by_detection_levels=True, by_cluster_expression=True)
vlm.normalize_by_total(min_perc_U=0.5)
vlm.adjust_totS_totU(normalize_total=True, fit_with_low_U=False, svr_C=1, svr_gamma=1e-04)
```


```
vlm.perform_PCA()
plt.plot(np.cumsum(vlm.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(vlm.pca.explained_variance_ratio_))>0.0055))[0][0]
vlm.pcs[:,1] *= -1 # flip for consistency with previous version
```


![png](kb_velocity_files/kb_velocity_57_0.png)



```
nn = NearestNeighbors(50)
nn.fit(vlm.pcs[:,:4])
knn_pca = nn.kneighbors_graph(mode='distance')
knn_pca = knn_pca.tocoo()
G = igraph.Graph(list(zip(knn_pca.row, knn_pca.col)), directed=False, edge_attrs={'weight': knn_pca.data})
VxCl = G.community_multilevel(return_levels=False, weights="weight")
labels = np.array(VxCl.membership)
```


```
pc_obj = principal_curve(vlm.pcs[:,:4], False)
pc_obj.arclength = np.max(pc_obj.arclength) - pc_obj.arclength
labels = np.argsort(np.argsort(aggregate_np(labels, pc_obj.arclength, func=np.median)))[labels]
```


```
manual_annotation = {str(i):[i] for i in labels}
annotation_dict = {v:k for k, values in manual_annotation.items() for v in values }
clusters = np.array([annotation_dict[i] for i in labels])
colors20 = np.vstack((plt.cm.tab20b(np.linspace(0., 1, 20))[::2], plt.cm.tab20c(np.linspace(0, 1, 20))[1::2]))
vlm.set_clusters(clusters, cluster_colors_dict={k:colors20[v[0] % 20,:] for k,v in manual_annotation.items()})
```


```
k = 550
vlm.knn_imputation(n_pca_dims=n_comps, k=k)
```


```
vlm.normalize_median()
vlm.fit_gammas(maxmin_perc=[2,95], limit_gamma=True)
```


```
vlm.normalize(which="imputed", size=False, log=True)
vlm.Pcs = np.array(vlm.pcs[:,:2], order="C")
```


```
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift()
vlm.extrapolate_cell_at_t(delta_t=1)
```


```
vlm.estimate_transition_prob(hidim="Sx_sz", embed="Pcs", transform="log", psc=1,
                             n_neighbors=150, knn_random=True, sampled_fraction=1)
```

    WARNING:root:Nans encountered in corrcoef and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.
    WARNING:root:Nans encountered in corrcoef_random and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.



```
vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=False)
```


```
vlm.calculate_grid_arrows(smooth=0.9, steps=(25, 25), n_neighbors=200)
```


```
plt.figure(None,(9,9))
vlm.plot_grid_arrows(scatter_kwargs_dict={"alpha":0.7, "lw":0.7, "edgecolor":"0.4", "s":70, "rasterized":True},
                     min_mass=2.9, angles='xy', scale_units='xy',
                     headaxislength=2.75, headlength=5, headwidth=4.8, quiver_scale=0.35, scale_type="absolute")
#plt.plot(pc_obj.projections[pc_obj.ixsort,0], pc_obj.projections[pc_obj.ixsort,1], c="w", lw=6, zorder=1000000)
#plt.plot(pc_obj.projections[pc_obj.ixsort,0], pc_obj.projections[pc_obj.ixsort,1], c="k", lw=3, zorder=2000000)
plt.gca().invert_xaxis()
plt.axis("off")
plt.axis("equal");
#plt.savefig("kallisto_velocity_forebrain_glut.pdf")
```

    WARNING:root:The arrow scale was set to be 'absolute' make sure you know how to properly interpret the plots



![png](kb_velocity_files/kb_velocity_68_1.png)

