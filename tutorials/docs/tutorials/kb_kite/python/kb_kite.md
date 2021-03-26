# Pre-processing and analysis of feature barcode single-cell RNA-seq data with KITE.

In this notebook, we will perform pre-processing and analysis of [10x Genomics pbmc_1k_protein_v3](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_protein_v3) feature barcoding dataset using the **Kallisto Indexing and Tag Extraction (KITE)** workflow, implemented with a wrapper called `kb`. It was developed by Kyung Hoi (Joseph) Min and A. Sina Booeshaghi.

In Feature Barcoding assays, cellular data are recorded as short DNA sequences using procedures adapted from single-cell RNA-seq. The **KITE** workflow generates a "Mismatch Map" containing the sequences of all Feature Barcodes used in the experiment as well as all of their single-base mismatches. The Mismatch Map is used to produce transcipt-to-gene (.t2g) and fasta (.fa) files to be used as inputs for kallisto. An index is made with kallisto index, then kallisto | bustools effectively searches the sequencing data for the sequences in the Mismatch Map.

## Pre-processing

### Download the data

__Note:__ We use the `-O` option for `wget` to rename the files to easily identify them.


```
%%time
!wget https://caltech.box.com/shared/static/asmj4nu90ydhsrk3pm7aaxu00cnnfige.txt -O checksums.txt
!wget https://caltech.box.com/shared/static/mp2vr3p6dztdyatuag8ir3cektmrztg8.gz -O pbmc_1k_protein_v3_antibody_S2_L001_R1_001.fastq.gz
!wget https://caltech.box.com/shared/static/f3payi1za7mn0jfai7vm10sy3yqwgpqh.gz -O pbmc_1k_protein_v3_antibody_S2_L001_R2_001.fastq.gz
!wget https://caltech.box.com/shared/static/e112bbczh9o1rl6gfin36bqp0ga7uvdy.gz -O pbmc_1k_protein_v3_antibody_S2_L002_R1_001.fastq.gz
!wget https://caltech.box.com/shared/static/3ve2axc8dr8v5nnrhmynrdgpqj6xg42k.gz -O pbmc_1k_protein_v3_antibody_S2_L002_R2_001.fastq.gz
```

    --2020-01-15 23:30:46--  https://caltech.box.com/shared/static/asmj4nu90ydhsrk3pm7aaxu00cnnfige.txt
    Resolving caltech.box.com (caltech.box.com)... 107.152.27.197, 107.152.26.197
    Connecting to caltech.box.com (caltech.box.com)|107.152.27.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/asmj4nu90ydhsrk3pm7aaxu00cnnfige.txt [following]
    --2020-01-15 23:30:47--  https://caltech.box.com/public/static/asmj4nu90ydhsrk3pm7aaxu00cnnfige.txt
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/asmj4nu90ydhsrk3pm7aaxu00cnnfige.txt [following]
    --2020-01-15 23:30:47--  https://caltech.app.box.com/public/static/asmj4nu90ydhsrk3pm7aaxu00cnnfige.txt
    Resolving caltech.app.box.com (caltech.app.box.com)... 107.152.25.199, 107.152.24.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|107.152.25.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!hB8g3gAf6IPFJZIeYelVLYaU0zgrjSKUuMGrh2XpaLbJfOaZbNvQJyU1WOhfEx3pQphmIKpk3tGVDAWoPtJEtOERAphLQ3pdMzYj4lOXNAbTVjZ-oA4CEHQCfXA5fFl8DH3LmHIDVMb0RwLhqtvPyTit97ekNENSEdrjW621kDkH0JFwovPzZ9iZPRNkhe-h3D7ccDgi7GMN5rMZkHbX5xy6VNoLlckTJE7atFuqVlH1M_71ryVd5kMgdZke9ohg34gaMzDGwPDyNyzfzV3x0tk5yHQTa1BlARCRamX-lWQyV272RNb6FKzlWI4UOrLB1C1Rq1IYAb1Ouaf6VcDvXWcrEB1DJSDRVULIHuqjTBeehAhYh2SStzO_h_rU2hZv7RG8LqUcWeMfQEW4Bc20wxj32BPyxh-EXTHvMdSbWWdLATy9zjYDOzQ1UNSodHKAd-LeFHnsdz0Sfs-3cO7-noW8kpuKfx8ZWeqYwphOnnLVPcDbJ-UHvvyRfxsShA3Xyc-QeSyjIyngkthDWnS-QGHDM6UsaMYjPsaeLeQJBn0zMRi1efjfNlNMILKeZOOSSA8yXonQ3HhtQn2ZLVVa30ESj-nA7GtK_XtrW7TFc8ga0k7JSx1JH4rIjkLhjUB53XDkO5Y3d3-xjwa7UwghZYaZMAay5F6NDoYMEZVo4iUrmT5DLwuVQOAfh1aqWkVmGAz5NgShbkfMI_nUAI5D6QBMPpHIIEB2QzJvEwVeDz5Us8JFUYhoCwv-BIJwoqxfeOv1qcNsb4Em2dF9zl-qnoDyz3taFfYyECyddu7UIDyRLVp938CfhMZkPAyAgw43QwGoIk76baaDzFJhWN3v3s7RHwECxC4DtxH41bPnR2IRoTGxxOdOyNvKnuHOwce29hdKS7MXsqltBnju99POYI-lqMNA4e9_ih3Mimvpf_ykTk2r_Cg3D1CnVUcXbRsLDz7LMmgLvLtDsqO4_bZrXNIdTnbh2g6Ad_t_zFELmsDPbdYrR-M2L9owyaKymj5FF63eAlW9KPTfBZRBCDgnW5hvn0Z4oqzjxnS54QLyu_rel-AThHArY0mYcEbuyJdmIpc9W4kCg2WNwL8Ewkzim37IVSZEtpZqi4T2ONHhKl2eoQFq8fdXZVcJ9bPHAHcSCk5TyTHwFvEovfQwIMCqhGxOGIOlI6uhtY7XYlarLYsjfsWDVS6F1iiwQDq3AnRE84RCIAcn6y_NrkFtMCDE-GsMbBy-cM0AvAy9T06Lf9RnYP4JzPBYmQ6FN8-C7d3ccy8VjtNBX44_dCkjdfKpDSJSngpnG4POHCDemJ8g51qEqvPMAKIu8Zo2OM9ehkr0YQd6tamuOMvBL0-UxMRCrD7H2Ox9WdiyXoRI7w../download [following]
    --2020-01-15 23:30:48--  https://public.boxcloud.com/d/1/b1!hB8g3gAf6IPFJZIeYelVLYaU0zgrjSKUuMGrh2XpaLbJfOaZbNvQJyU1WOhfEx3pQphmIKpk3tGVDAWoPtJEtOERAphLQ3pdMzYj4lOXNAbTVjZ-oA4CEHQCfXA5fFl8DH3LmHIDVMb0RwLhqtvPyTit97ekNENSEdrjW621kDkH0JFwovPzZ9iZPRNkhe-h3D7ccDgi7GMN5rMZkHbX5xy6VNoLlckTJE7atFuqVlH1M_71ryVd5kMgdZke9ohg34gaMzDGwPDyNyzfzV3x0tk5yHQTa1BlARCRamX-lWQyV272RNb6FKzlWI4UOrLB1C1Rq1IYAb1Ouaf6VcDvXWcrEB1DJSDRVULIHuqjTBeehAhYh2SStzO_h_rU2hZv7RG8LqUcWeMfQEW4Bc20wxj32BPyxh-EXTHvMdSbWWdLATy9zjYDOzQ1UNSodHKAd-LeFHnsdz0Sfs-3cO7-noW8kpuKfx8ZWeqYwphOnnLVPcDbJ-UHvvyRfxsShA3Xyc-QeSyjIyngkthDWnS-QGHDM6UsaMYjPsaeLeQJBn0zMRi1efjfNlNMILKeZOOSSA8yXonQ3HhtQn2ZLVVa30ESj-nA7GtK_XtrW7TFc8ga0k7JSx1JH4rIjkLhjUB53XDkO5Y3d3-xjwa7UwghZYaZMAay5F6NDoYMEZVo4iUrmT5DLwuVQOAfh1aqWkVmGAz5NgShbkfMI_nUAI5D6QBMPpHIIEB2QzJvEwVeDz5Us8JFUYhoCwv-BIJwoqxfeOv1qcNsb4Em2dF9zl-qnoDyz3taFfYyECyddu7UIDyRLVp938CfhMZkPAyAgw43QwGoIk76baaDzFJhWN3v3s7RHwECxC4DtxH41bPnR2IRoTGxxOdOyNvKnuHOwce29hdKS7MXsqltBnju99POYI-lqMNA4e9_ih3Mimvpf_ykTk2r_Cg3D1CnVUcXbRsLDz7LMmgLvLtDsqO4_bZrXNIdTnbh2g6Ad_t_zFELmsDPbdYrR-M2L9owyaKymj5FF63eAlW9KPTfBZRBCDgnW5hvn0Z4oqzjxnS54QLyu_rel-AThHArY0mYcEbuyJdmIpc9W4kCg2WNwL8Ewkzim37IVSZEtpZqi4T2ONHhKl2eoQFq8fdXZVcJ9bPHAHcSCk5TyTHwFvEovfQwIMCqhGxOGIOlI6uhtY7XYlarLYsjfsWDVS6F1iiwQDq3AnRE84RCIAcn6y_NrkFtMCDE-GsMbBy-cM0AvAy9T06Lf9RnYP4JzPBYmQ6FN8-C7d3ccy8VjtNBX44_dCkjdfKpDSJSngpnG4POHCDemJ8g51qEqvPMAKIu8Zo2OM9ehkr0YQd6tamuOMvBL0-UxMRCrD7H2Ox9WdiyXoRI7w../download
    Resolving public.boxcloud.com (public.boxcloud.com)... 107.152.25.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|107.152.25.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 510 [text/plain]
    Saving to: ‘checksums.txt’
    
    checksums.txt       100%[===================>]     510  --.-KB/s    in 0s      
    
    2020-01-15 23:30:48 (48.5 MB/s) - ‘checksums.txt’ saved [510/510]
    
    --2020-01-15 23:30:49--  https://caltech.box.com/shared/static/mp2vr3p6dztdyatuag8ir3cektmrztg8.gz
    Resolving caltech.box.com (caltech.box.com)... 107.152.27.197, 107.152.26.197
    Connecting to caltech.box.com (caltech.box.com)|107.152.27.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/mp2vr3p6dztdyatuag8ir3cektmrztg8.gz [following]
    --2020-01-15 23:30:49--  https://caltech.box.com/public/static/mp2vr3p6dztdyatuag8ir3cektmrztg8.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/mp2vr3p6dztdyatuag8ir3cektmrztg8.gz [following]
    --2020-01-15 23:30:49--  https://caltech.app.box.com/public/static/mp2vr3p6dztdyatuag8ir3cektmrztg8.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 107.152.25.199, 107.152.24.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|107.152.25.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!wC_RYxsJHveXuH8OpwfnHoJmeUtNtr6DwJGYhcQl1woKrwGO7wLfQcu2Sj7_r_5PP0EMBhI04ARBvSnoQ3KisUsiTrdNdBRS4JI9Kb83gDBI2lV0Z-1x7foync7tehoEPKavSVhvgakGHbLGnKyjnFEDPbbZT-xAZKzvK7kT6rX-0zc7wWJNzO72SGvsTgPLe1Z18Oyd4Nr4NeWM2Cjz7xfnDeINRQvhoCOe_PfFaluBu3hZsWjQSTAhr20w5nS8udP9brb02MsB3xkmgxls55acTrJyhlH2z7pzjrWKzxWmOWi3wDyelw66oDRtXJ5I-X5Gp73yC6UTBjrqX0-Q-BMdQ8Fb7Y3FWXgtKiNiJqYjMA-hwJDXl6edc_DGNq3amal01GrOj7mY8NWK8a-iMJt66jC7pP72i4KQOSOyhC6mr5pThvYmrbfBCT3l5IL2ynQoGRnAuO9PqnRSpHv2JTx95TJVllgy42byPkpu9PQv1tFyrS8UfpWHVR4aEWC6SmdnaZ7mqSBXQNQMs1ra8qY8PIRVvUR13uo57qlLNbB6jBQK3QznwFPfWqst04YpI1VcsIRlU72JPd1ukDvpC3F51OyMIsQNeaQsbIUDVt_k7kOkg5pykLWdwMfPqvSWHx7uCfUPmZZmoiLomtLFNuwxXMOKfN2AiaggGMTfKs4vljsiFfjFjPpNoVaunEKiZ4FTTrMip4I7u2bvzzz-K37e8P1OnZWhb7-pEyoNWzS1uyOh7HtCQOT11bEMKde1rQr3wQt_vRauhUh7xVQ6wZ-ao1OGUs3G7jgb93khZJSKB2zpedD-YmRIvla1cTYV6SLDkXjGqDJjpPhQ5YaiMVR0x_5UZJ9mu1h5UtQucKO-f35JQuXrNbG6pSXlypsYrPokaDR6jzglS3Px7kDDZSD1Im5O0kDf7D9LTDxF_wzrVFJdOIUfAJXokRMwctKaotTFzDWP9Dy1lBm8lPR82igNge-a-uEqF8r1jEb5tVknD-z6z0VGy2Bn4oa941pZEdx1wBwOX6mjZktDyCHuCmkkZ2fRPSu8Qc90Jceth8RiQVZHX5Q6sefKH3nUusgrs0hfyevGnmQplyatNnwfvrSstco7bjkQwo7-PxfHk4y38SbWIOE_cAocNnDDlKOMzlrs3cKXx8-RW06g3VWpJgQUXJq6f2X61Af05NEM6PNh_z2kKB_M8XueOjvIYflWHkaOIyj8XsPjdakI3_4uyIoZgxo_B1NBcTwgHIJxdW9zUmG7fpBQvv6_Z8ueHHIMPdiSOT6wD-tofkBUnprGWlrBiy2H9fBgnEzx7EwjUfaNXQ5CUSNH0uu2nyki6YrGkZ7FCYOTCyxOooqdF932nyEz2IXf_VDADzg8QIDX4mV5j6k2tGY3QaIfHbioteyEDx7zJcAlzB_Y9ddZHsjNrAxwlXTBPHYlRwQm4jAYmu4ZMutTBKo./download [following]
    --2020-01-15 23:30:50--  https://public.boxcloud.com/d/1/b1!wC_RYxsJHveXuH8OpwfnHoJmeUtNtr6DwJGYhcQl1woKrwGO7wLfQcu2Sj7_r_5PP0EMBhI04ARBvSnoQ3KisUsiTrdNdBRS4JI9Kb83gDBI2lV0Z-1x7foync7tehoEPKavSVhvgakGHbLGnKyjnFEDPbbZT-xAZKzvK7kT6rX-0zc7wWJNzO72SGvsTgPLe1Z18Oyd4Nr4NeWM2Cjz7xfnDeINRQvhoCOe_PfFaluBu3hZsWjQSTAhr20w5nS8udP9brb02MsB3xkmgxls55acTrJyhlH2z7pzjrWKzxWmOWi3wDyelw66oDRtXJ5I-X5Gp73yC6UTBjrqX0-Q-BMdQ8Fb7Y3FWXgtKiNiJqYjMA-hwJDXl6edc_DGNq3amal01GrOj7mY8NWK8a-iMJt66jC7pP72i4KQOSOyhC6mr5pThvYmrbfBCT3l5IL2ynQoGRnAuO9PqnRSpHv2JTx95TJVllgy42byPkpu9PQv1tFyrS8UfpWHVR4aEWC6SmdnaZ7mqSBXQNQMs1ra8qY8PIRVvUR13uo57qlLNbB6jBQK3QznwFPfWqst04YpI1VcsIRlU72JPd1ukDvpC3F51OyMIsQNeaQsbIUDVt_k7kOkg5pykLWdwMfPqvSWHx7uCfUPmZZmoiLomtLFNuwxXMOKfN2AiaggGMTfKs4vljsiFfjFjPpNoVaunEKiZ4FTTrMip4I7u2bvzzz-K37e8P1OnZWhb7-pEyoNWzS1uyOh7HtCQOT11bEMKde1rQr3wQt_vRauhUh7xVQ6wZ-ao1OGUs3G7jgb93khZJSKB2zpedD-YmRIvla1cTYV6SLDkXjGqDJjpPhQ5YaiMVR0x_5UZJ9mu1h5UtQucKO-f35JQuXrNbG6pSXlypsYrPokaDR6jzglS3Px7kDDZSD1Im5O0kDf7D9LTDxF_wzrVFJdOIUfAJXokRMwctKaotTFzDWP9Dy1lBm8lPR82igNge-a-uEqF8r1jEb5tVknD-z6z0VGy2Bn4oa941pZEdx1wBwOX6mjZktDyCHuCmkkZ2fRPSu8Qc90Jceth8RiQVZHX5Q6sefKH3nUusgrs0hfyevGnmQplyatNnwfvrSstco7bjkQwo7-PxfHk4y38SbWIOE_cAocNnDDlKOMzlrs3cKXx8-RW06g3VWpJgQUXJq6f2X61Af05NEM6PNh_z2kKB_M8XueOjvIYflWHkaOIyj8XsPjdakI3_4uyIoZgxo_B1NBcTwgHIJxdW9zUmG7fpBQvv6_Z8ueHHIMPdiSOT6wD-tofkBUnprGWlrBiy2H9fBgnEzx7EwjUfaNXQ5CUSNH0uu2nyki6YrGkZ7FCYOTCyxOooqdF932nyEz2IXf_VDADzg8QIDX4mV5j6k2tGY3QaIfHbioteyEDx7zJcAlzB_Y9ddZHsjNrAxwlXTBPHYlRwQm4jAYmu4ZMutTBKo./download
    Resolving public.boxcloud.com (public.boxcloud.com)... 107.152.25.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|107.152.25.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 137613323 (131M) [application/octet-stream]
    Saving to: ‘pbmc_1k_protein_v3_antibody_S2_L001_R1_001.fastq.gz’
    
    pbmc_1k_protein_v3_ 100%[===================>] 131.24M  20.7MB/s    in 6.6s    
    
    2020-01-15 23:30:57 (19.8 MB/s) - ‘pbmc_1k_protein_v3_antibody_S2_L001_R1_001.fastq.gz’ saved [137613323/137613323]
    
    --2020-01-15 23:30:57--  https://caltech.box.com/shared/static/f3payi1za7mn0jfai7vm10sy3yqwgpqh.gz
    Resolving caltech.box.com (caltech.box.com)... 107.152.27.197, 107.152.26.197
    Connecting to caltech.box.com (caltech.box.com)|107.152.27.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/f3payi1za7mn0jfai7vm10sy3yqwgpqh.gz [following]
    --2020-01-15 23:30:58--  https://caltech.box.com/public/static/f3payi1za7mn0jfai7vm10sy3yqwgpqh.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/f3payi1za7mn0jfai7vm10sy3yqwgpqh.gz [following]
    --2020-01-15 23:30:58--  https://caltech.app.box.com/public/static/f3payi1za7mn0jfai7vm10sy3yqwgpqh.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 107.152.25.199, 107.152.24.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|107.152.25.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!fBmz4pkQcuY_wkH_de3fFELUuflKEu6dQTmRu68vex9y0w8J-9-Sxl2WVtUr2YTIfjcmcrLYtIGqELWqsd-acgDGkR6_UZ0ZXCeUhz8X2yRO3qZMA9NHadT2BuyJxRP34tsakvziyemt-JCQGL32VjRl3yJhcKvKBIp1IMZT5LytUuUhNAdkc_uyyQ_MweXWXKFxsXaKlk8vtDs4j9K2MQr_q8LqGF8Fc0hGm222OmB_sz41sjXMDFAC4MrMq2RrtMm7pIt0HNLBq2yL554O2IGLIVidM_AdJho5yQZKDVG40BjaCb8GPvtuRedw_gWJthGO82PEJNAlTfXO2lwDyX2nK5A78gwu1s8BZrsO_0bNYT-bsefyKcf3QUtrcl8oQAsoPbjyeOcim5L1u1JcQ67Bm3EqXGn4A6k14pDVdemUtoNLmwsAtGBIpFnYLa7l4L-URz7uMQw8KYTg6PFj3TU96HMWUTOOUlxBtBysW4SsiJLlCERSM9VAnHkPxFE6QJrrG5u8IUNQkjTaMZ9Bf_rANq_C0cIG3S9R1VOaypdF3YFUlmufKFQStUHjZ3nfzRwChcyDinP0YOqCNvshNJ1lzxTfWsXE2OdR05nhwfwF64sYGhZRxO9dyrw9USRtxBAl38JwNmcbeZ8_xXVelJN3jn7Ograq1bFogoZUTZ5rhgcEF3Er5AIrZRmWoEDPuhpTHEzljNRQQZQvklf9NaeJu2x4_8xQz_seyiw832AZm-i-XoWSYHPEyPHVKyqUrRZiJEcNlZher94wpg03fjhJsZtsBmid9KmR99xGuP27_CcKWKcUzJFhPaJ2FO86349g50POwVlUVUBkdmCZGn1-boGArraPu4Uojg9yE7kdHOKxC7nLN4DmiTVg-sqZWcV26asJwthTnwG-hei-FQtlZwK4Sqpx5IVHnQDh5hjxhw97UZhjulsa2ub5PgXMaULyHZ1e5-Cc2ZKTgvMvfesP-rmebI0FSc3rDA58TUkUAXbwPnrFBknLhKItYF-TBuXm0i0E5Mc5DdgI7R8loqSZK9c5KreSS61UhXxFVmBVw5zQ_vWqtz7JSwQDdX4LhROLotFcog8rtDPN6cHzFG8XRen4KGi45AkMiMQG8Qd-9H8eTYGQevynMmky_ikO-pMGwjkINWL4LxrnOEqrLXN9zVaSvoh0DkLHhGQtdvP4aFS0fzhEU7HYgjP17ONhK5MRzQa_USiryd5_49xNLTG7qr3BlX6QC6-6Lt4mGDSMETvpuVBRkotupHYCRQUDGEGqNFCTOvjOmvEf3AG_xFC9NQaB4tgTaQ7ct5-muXR2quScvPwT62N1jTvCKBNmm-fC5xG3ioykZaf9JHtJ8hWpqZtiGbQBz9WdnTlOmpHgpnat22EEgKjhWYJuWB5N4NZYisDK-kQT3cofhDnDelGd5exguZWjyfV1vnT5-Ws0OTkA4w../download [following]
    --2020-01-15 23:30:58--  https://public.boxcloud.com/d/1/b1!fBmz4pkQcuY_wkH_de3fFELUuflKEu6dQTmRu68vex9y0w8J-9-Sxl2WVtUr2YTIfjcmcrLYtIGqELWqsd-acgDGkR6_UZ0ZXCeUhz8X2yRO3qZMA9NHadT2BuyJxRP34tsakvziyemt-JCQGL32VjRl3yJhcKvKBIp1IMZT5LytUuUhNAdkc_uyyQ_MweXWXKFxsXaKlk8vtDs4j9K2MQr_q8LqGF8Fc0hGm222OmB_sz41sjXMDFAC4MrMq2RrtMm7pIt0HNLBq2yL554O2IGLIVidM_AdJho5yQZKDVG40BjaCb8GPvtuRedw_gWJthGO82PEJNAlTfXO2lwDyX2nK5A78gwu1s8BZrsO_0bNYT-bsefyKcf3QUtrcl8oQAsoPbjyeOcim5L1u1JcQ67Bm3EqXGn4A6k14pDVdemUtoNLmwsAtGBIpFnYLa7l4L-URz7uMQw8KYTg6PFj3TU96HMWUTOOUlxBtBysW4SsiJLlCERSM9VAnHkPxFE6QJrrG5u8IUNQkjTaMZ9Bf_rANq_C0cIG3S9R1VOaypdF3YFUlmufKFQStUHjZ3nfzRwChcyDinP0YOqCNvshNJ1lzxTfWsXE2OdR05nhwfwF64sYGhZRxO9dyrw9USRtxBAl38JwNmcbeZ8_xXVelJN3jn7Ograq1bFogoZUTZ5rhgcEF3Er5AIrZRmWoEDPuhpTHEzljNRQQZQvklf9NaeJu2x4_8xQz_seyiw832AZm-i-XoWSYHPEyPHVKyqUrRZiJEcNlZher94wpg03fjhJsZtsBmid9KmR99xGuP27_CcKWKcUzJFhPaJ2FO86349g50POwVlUVUBkdmCZGn1-boGArraPu4Uojg9yE7kdHOKxC7nLN4DmiTVg-sqZWcV26asJwthTnwG-hei-FQtlZwK4Sqpx5IVHnQDh5hjxhw97UZhjulsa2ub5PgXMaULyHZ1e5-Cc2ZKTgvMvfesP-rmebI0FSc3rDA58TUkUAXbwPnrFBknLhKItYF-TBuXm0i0E5Mc5DdgI7R8loqSZK9c5KreSS61UhXxFVmBVw5zQ_vWqtz7JSwQDdX4LhROLotFcog8rtDPN6cHzFG8XRen4KGi45AkMiMQG8Qd-9H8eTYGQevynMmky_ikO-pMGwjkINWL4LxrnOEqrLXN9zVaSvoh0DkLHhGQtdvP4aFS0fzhEU7HYgjP17ONhK5MRzQa_USiryd5_49xNLTG7qr3BlX6QC6-6Lt4mGDSMETvpuVBRkotupHYCRQUDGEGqNFCTOvjOmvEf3AG_xFC9NQaB4tgTaQ7ct5-muXR2quScvPwT62N1jTvCKBNmm-fC5xG3ioykZaf9JHtJ8hWpqZtiGbQBz9WdnTlOmpHgpnat22EEgKjhWYJuWB5N4NZYisDK-kQT3cofhDnDelGd5exguZWjyfV1vnT5-Ws0OTkA4w../download
    Resolving public.boxcloud.com (public.boxcloud.com)... 107.152.25.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|107.152.25.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 268437558 (256M) [application/octet-stream]
    Saving to: ‘pbmc_1k_protein_v3_antibody_S2_L001_R2_001.fastq.gz’
    
    pbmc_1k_protein_v3_ 100%[===================>] 256.00M  25.6MB/s    in 10s     
    
    2020-01-15 23:31:09 (24.9 MB/s) - ‘pbmc_1k_protein_v3_antibody_S2_L001_R2_001.fastq.gz’ saved [268437558/268437558]
    
    --2020-01-15 23:31:10--  https://caltech.box.com/shared/static/e112bbczh9o1rl6gfin36bqp0ga7uvdy.gz
    Resolving caltech.box.com (caltech.box.com)... 107.152.27.197, 107.152.26.197
    Connecting to caltech.box.com (caltech.box.com)|107.152.27.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/e112bbczh9o1rl6gfin36bqp0ga7uvdy.gz [following]
    --2020-01-15 23:31:10--  https://caltech.box.com/public/static/e112bbczh9o1rl6gfin36bqp0ga7uvdy.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/e112bbczh9o1rl6gfin36bqp0ga7uvdy.gz [following]
    --2020-01-15 23:31:10--  https://caltech.app.box.com/public/static/e112bbczh9o1rl6gfin36bqp0ga7uvdy.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 107.152.25.199, 107.152.24.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|107.152.25.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!M59nnZqPulpKuBWtry8EHWYYUiaaOsIKjWqHd0IQXy7aruqzgn-CGHSOl8B12Db_a_grwP4zcLOEP5IE-owI-xeJPjEer1_anWScWNU60iE3SqaoBOWgMKhcCex4nDj8VgWI2dAOU6qV9ZWB5J8QuhoVw2ae9ZPKu-jwGYup6rfCxYqPDkCAf0lnz7XRoSm4ggYwxYgy-HpRcfJ2GgDaYITwuqboLejkZr6NleBGygj-PXdgtMDkv9TZdUh2JfkkpIFMeQT6szx9rpC8B9fBB88XAGHKQonNNUSco6gT--kfU6oqwY27wino5DXNW3xZKU93FUqlIFGPb8OlulzX26cee9rN7qtYeYaVd5uX6x4yuh-nDiPYwkHv9rA1Wjz9rXUnCNpEpnl1-YhzKSyST_6jlrWbPDNoSBr8UC1Ox_drhYcXv5pUykgazaJ543TaHTi1pfBWKJC9GQHjLxtgcwxfwg94dyPrRpbzmaMczYrlNwkWXf7h89w0YFEUfU7oFyLXsDZhgdpLiNTNJgzjvN-tjpJeHQtD258NrKfKs5C4sTZZ-9QlJLREOQwCtBXs5He8HaLcMeUnTyl9R5rVBXTxbZJXqcN2GRTIR4Hr1WEbDy14JcX7PsxlkOr4_scjHFKhiT9SzqVZEOumHUHQrUMfeqYxpnje3kcFe-HW-tK4dAW4_De8OLn81G2Iluhk31-DcITZSrLcNcTWExkFnLxTgUaS7HUBlOqUkuTK20QUuKWzGsBSdAOOOK2gXcMetPNMPdN_f_TF_TaQh8-kZRTFSlwklcEYXSBk9tSvya6g8I73fx0HY2CLYV_1wm-I50gbwkJ-B73sEfQVU2LHb8ZiWYOCIf7Mva19bJMwli9_EIlC9GX52DAlp_iGiBuP3vpmVuUrgl7mKd1DfgD1rZVIpBbZlGYRc_52tKn9LarF6Oo0K3e66bWvH92xoORwV7sgOAw-vg5cVZcYJrzbOnJS5-Uo3cVz3ZpDw4DaZRVngLSFMPYEhkYGKrBFxJF1TuU6DP0-_0WbG6Gj_cG8GCMttfomu0V7d3E5H1PeVsaTB2ZHzACp0fyaRpK5mjseGicTF9n8l_lFRTuGJ3xXC10n0p3rHbRG2dV4VFWAkvot-to2PpGsmfG_ZOtj_bEBM5SK_JfO5Ku7gJcFlMFMRP7rZGTaDG9xevkwMythSR-L_nrNoWwlkHzdJQJI-p7um1S3U_RlTC7cu2ioZXFWD8AQIM2e2EXj3u0IvAmoIOu6yXF8QFGjzQQAZgS4PJZfG5CdIJzPazkLp6eQR30TfXjrpwmLPaQgaZwLkFYNjUqii6iFQd6tOmEzKg9OniXJev0_oia5Rc0FU1zlKU-Mfk6W6muyrzjeg31M_i33Q_zx-g22-ZoaSZLdiacFXo2hY83YtSOTKZUANL8oDX8f2EE-79gcsMRTJdag4bTCc4k9K3sNPQ../download [following]
    --2020-01-15 23:31:11--  https://public.boxcloud.com/d/1/b1!M59nnZqPulpKuBWtry8EHWYYUiaaOsIKjWqHd0IQXy7aruqzgn-CGHSOl8B12Db_a_grwP4zcLOEP5IE-owI-xeJPjEer1_anWScWNU60iE3SqaoBOWgMKhcCex4nDj8VgWI2dAOU6qV9ZWB5J8QuhoVw2ae9ZPKu-jwGYup6rfCxYqPDkCAf0lnz7XRoSm4ggYwxYgy-HpRcfJ2GgDaYITwuqboLejkZr6NleBGygj-PXdgtMDkv9TZdUh2JfkkpIFMeQT6szx9rpC8B9fBB88XAGHKQonNNUSco6gT--kfU6oqwY27wino5DXNW3xZKU93FUqlIFGPb8OlulzX26cee9rN7qtYeYaVd5uX6x4yuh-nDiPYwkHv9rA1Wjz9rXUnCNpEpnl1-YhzKSyST_6jlrWbPDNoSBr8UC1Ox_drhYcXv5pUykgazaJ543TaHTi1pfBWKJC9GQHjLxtgcwxfwg94dyPrRpbzmaMczYrlNwkWXf7h89w0YFEUfU7oFyLXsDZhgdpLiNTNJgzjvN-tjpJeHQtD258NrKfKs5C4sTZZ-9QlJLREOQwCtBXs5He8HaLcMeUnTyl9R5rVBXTxbZJXqcN2GRTIR4Hr1WEbDy14JcX7PsxlkOr4_scjHFKhiT9SzqVZEOumHUHQrUMfeqYxpnje3kcFe-HW-tK4dAW4_De8OLn81G2Iluhk31-DcITZSrLcNcTWExkFnLxTgUaS7HUBlOqUkuTK20QUuKWzGsBSdAOOOK2gXcMetPNMPdN_f_TF_TaQh8-kZRTFSlwklcEYXSBk9tSvya6g8I73fx0HY2CLYV_1wm-I50gbwkJ-B73sEfQVU2LHb8ZiWYOCIf7Mva19bJMwli9_EIlC9GX52DAlp_iGiBuP3vpmVuUrgl7mKd1DfgD1rZVIpBbZlGYRc_52tKn9LarF6Oo0K3e66bWvH92xoORwV7sgOAw-vg5cVZcYJrzbOnJS5-Uo3cVz3ZpDw4DaZRVngLSFMPYEhkYGKrBFxJF1TuU6DP0-_0WbG6Gj_cG8GCMttfomu0V7d3E5H1PeVsaTB2ZHzACp0fyaRpK5mjseGicTF9n8l_lFRTuGJ3xXC10n0p3rHbRG2dV4VFWAkvot-to2PpGsmfG_ZOtj_bEBM5SK_JfO5Ku7gJcFlMFMRP7rZGTaDG9xevkwMythSR-L_nrNoWwlkHzdJQJI-p7um1S3U_RlTC7cu2ioZXFWD8AQIM2e2EXj3u0IvAmoIOu6yXF8QFGjzQQAZgS4PJZfG5CdIJzPazkLp6eQR30TfXjrpwmLPaQgaZwLkFYNjUqii6iFQd6tOmEzKg9OniXJev0_oia5Rc0FU1zlKU-Mfk6W6muyrzjeg31M_i33Q_zx-g22-ZoaSZLdiacFXo2hY83YtSOTKZUANL8oDX8f2EE-79gcsMRTJdag4bTCc4k9K3sNPQ../download
    Resolving public.boxcloud.com (public.boxcloud.com)... 107.152.25.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|107.152.25.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 133963077 (128M) [application/octet-stream]
    Saving to: ‘pbmc_1k_protein_v3_antibody_S2_L002_R1_001.fastq.gz’
    
    pbmc_1k_protein_v3_ 100%[===================>] 127.76M  19.8MB/s    in 6.7s    
    
    2020-01-15 23:31:18 (19.1 MB/s) - ‘pbmc_1k_protein_v3_antibody_S2_L002_R1_001.fastq.gz’ saved [133963077/133963077]
    
    --2020-01-15 23:31:19--  https://caltech.box.com/shared/static/3ve2axc8dr8v5nnrhmynrdgpqj6xg42k.gz
    Resolving caltech.box.com (caltech.box.com)... 107.152.27.197, 107.152.26.197
    Connecting to caltech.box.com (caltech.box.com)|107.152.27.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/3ve2axc8dr8v5nnrhmynrdgpqj6xg42k.gz [following]
    --2020-01-15 23:31:19--  https://caltech.box.com/public/static/3ve2axc8dr8v5nnrhmynrdgpqj6xg42k.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/3ve2axc8dr8v5nnrhmynrdgpqj6xg42k.gz [following]
    --2020-01-15 23:31:19--  https://caltech.app.box.com/public/static/3ve2axc8dr8v5nnrhmynrdgpqj6xg42k.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 107.152.25.199, 107.152.24.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|107.152.25.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!1o59St-Vp1K_5H6iqjyqoRtxUre8NDxVi6vrJEXQ5gByoCHBMqaoteEOl71swIhN-m8I5iJeFMi07M4qsnm7ktXV2qDfxrUBFwMZXCps80uDraGaNkZeuy1r7BteFZ47vPPe1PRyyxaJT1uM2KKIUFKegjO5KWmqkAR7Og7X0gDDPoo9QLzqTACrlD_2TFO8TwVNRjhfjMlgrjj1CC8n3QZ-_wCr2VqGnCkhLQSsZ1bAUKbcV_O4-X9YcXTFl1aEc8QHWE8cCoQ-WhqUSJ2P-tx1eM6BLvg2pXkjn_eYL0Di5-_RnxX4ilKsomSrsF7lYYvTuBVlmtvZiAiDp7VhdOX40jJTAMAKE6BgFlondsID2wXlNuM_TKqWzBquwv5pRbm9Lw7YBcZoDjAz5fSAnhORsiwcm-pkkBf9Qwr0r5zHLz_dD3xGdFxEpjssOL5LYWng-QW7WZMWqT3IpQUtgL3YMtFEfBSc5V3xU9V7nklR4kM6OoASYS7Rjmu_JuTnBtxpdW1U3EtU0GduqsbuYVpteV2EPZjvNcOgbBN5OBGzkQYAaeKA9s4iAf6FWPH69xS865QaXulpneiDxm4f7UK0UoVytXS4sEd4eJAlLT4YUGInUI56cdrhUhtIKxEdOTCbs7eOOC-SWtvRQTLvzu0dY_JwPmU3IJhzKThRqs8Oq2VJjYNRaoljHVOR-8SV7dqkkI-4WJvFThndsygpweHpqpUVfCNaeikbUcwcAy0OM4ZewR80IfGfnGcyrwzWkEZm9LLqr7pU0kQ2flPjngTBzQ__C5OsX7Z_YdDaLOetDhqh0JJI3KkObvVs8XzXtdcAXudqfH3V6PsfizH2C85__YmpLu88FTqTdaJzfwvOr8dNt7dbx1mA-gkbUxJK1Hpzt2beu_DL23UjZA0n8dA_RfVVbDETwWRqOEFtvQPgDS7TLIAEzMu1MoLSJE4uuM-O2kSff5rGbx9dJEEjx4qSZK4xLouiU5wy1a2jlbDE6bDx8tj5vlZrqll0x8NAIwJMGa4qBScbI5IKM3k-Zi0b5kuxd5vWOlLr8diIUOwh06GSNHMQcZCL2abT69eEf0OrJerVaPxPbVPpqUQmY96NGA7TcAtpBfxwWRZH0rFrJM1_6gyiMSiyMEXdr3-Zf-7foZeeWIT2pBdTg_xUyAj_9TthjNf0e3LUK8KqpAhuKVqAXYHdT2plLaW9NjtfytbpDXv8CoJH_HLx2NNi-4zhz_-PC6hIDgVWJevwsLHtScmLghgNCwiRMDyBgCeXBdFubPlMWES5Zhb_clsQfOBQShMvkuom1nDBREdxjiHXHHazDN2KIsFBG1iX3KmYqH3s50Kpgm1WeD3OIrT6LQ4m3boVu4Cg7goyXdBHuRnbxLDaUt1EYDw7Wi4_yYyTQmsIHdfJ-jR1S8WXtkttxw03NxGoh2LWxSL7Dt9apr2RjlRcers./download [following]
    --2020-01-15 23:31:20--  https://public.boxcloud.com/d/1/b1!1o59St-Vp1K_5H6iqjyqoRtxUre8NDxVi6vrJEXQ5gByoCHBMqaoteEOl71swIhN-m8I5iJeFMi07M4qsnm7ktXV2qDfxrUBFwMZXCps80uDraGaNkZeuy1r7BteFZ47vPPe1PRyyxaJT1uM2KKIUFKegjO5KWmqkAR7Og7X0gDDPoo9QLzqTACrlD_2TFO8TwVNRjhfjMlgrjj1CC8n3QZ-_wCr2VqGnCkhLQSsZ1bAUKbcV_O4-X9YcXTFl1aEc8QHWE8cCoQ-WhqUSJ2P-tx1eM6BLvg2pXkjn_eYL0Di5-_RnxX4ilKsomSrsF7lYYvTuBVlmtvZiAiDp7VhdOX40jJTAMAKE6BgFlondsID2wXlNuM_TKqWzBquwv5pRbm9Lw7YBcZoDjAz5fSAnhORsiwcm-pkkBf9Qwr0r5zHLz_dD3xGdFxEpjssOL5LYWng-QW7WZMWqT3IpQUtgL3YMtFEfBSc5V3xU9V7nklR4kM6OoASYS7Rjmu_JuTnBtxpdW1U3EtU0GduqsbuYVpteV2EPZjvNcOgbBN5OBGzkQYAaeKA9s4iAf6FWPH69xS865QaXulpneiDxm4f7UK0UoVytXS4sEd4eJAlLT4YUGInUI56cdrhUhtIKxEdOTCbs7eOOC-SWtvRQTLvzu0dY_JwPmU3IJhzKThRqs8Oq2VJjYNRaoljHVOR-8SV7dqkkI-4WJvFThndsygpweHpqpUVfCNaeikbUcwcAy0OM4ZewR80IfGfnGcyrwzWkEZm9LLqr7pU0kQ2flPjngTBzQ__C5OsX7Z_YdDaLOetDhqh0JJI3KkObvVs8XzXtdcAXudqfH3V6PsfizH2C85__YmpLu88FTqTdaJzfwvOr8dNt7dbx1mA-gkbUxJK1Hpzt2beu_DL23UjZA0n8dA_RfVVbDETwWRqOEFtvQPgDS7TLIAEzMu1MoLSJE4uuM-O2kSff5rGbx9dJEEjx4qSZK4xLouiU5wy1a2jlbDE6bDx8tj5vlZrqll0x8NAIwJMGa4qBScbI5IKM3k-Zi0b5kuxd5vWOlLr8diIUOwh06GSNHMQcZCL2abT69eEf0OrJerVaPxPbVPpqUQmY96NGA7TcAtpBfxwWRZH0rFrJM1_6gyiMSiyMEXdr3-Zf-7foZeeWIT2pBdTg_xUyAj_9TthjNf0e3LUK8KqpAhuKVqAXYHdT2plLaW9NjtfytbpDXv8CoJH_HLx2NNi-4zhz_-PC6hIDgVWJevwsLHtScmLghgNCwiRMDyBgCeXBdFubPlMWES5Zhb_clsQfOBQShMvkuom1nDBREdxjiHXHHazDN2KIsFBG1iX3KmYqH3s50Kpgm1WeD3OIrT6LQ4m3boVu4Cg7goyXdBHuRnbxLDaUt1EYDw7Wi4_yYyTQmsIHdfJ-jR1S8WXtkttxw03NxGoh2LWxSL7Dt9apr2RjlRcers./download
    Resolving public.boxcloud.com (public.boxcloud.com)... 107.152.25.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|107.152.25.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 262167194 (250M) [application/octet-stream]
    Saving to: ‘pbmc_1k_protein_v3_antibody_S2_L002_R2_001.fastq.gz’
    
    pbmc_1k_protein_v3_ 100%[===================>] 250.02M  23.2MB/s    in 11s     
    
    2020-01-15 23:31:31 (22.9 MB/s) - ‘pbmc_1k_protein_v3_antibody_S2_L002_R2_001.fastq.gz’ saved [262167194/262167194]
    
    CPU times: user 455 ms, sys: 85.5 ms, total: 540 ms
    Wall time: 45.4 s


Then, we verify the integrity of the files we downloaded to make sure they were not corrupted during the download.


```
!md5sum -c checksums.txt --ignore-missing
```

    pbmc_1k_protein_v3_antibody_S2_L001_R1_001.fastq.gz: OK
    pbmc_1k_protein_v3_antibody_S2_L001_R2_001.fastq.gz: OK
    pbmc_1k_protein_v3_antibody_S2_L002_R1_001.fastq.gz: OK
    pbmc_1k_protein_v3_antibody_S2_L002_R2_001.fastq.gz: OK


### Install `kb`

Install `kb` for running the kallisto|bustools workflow.


```
!pip install git+https://github.com/pachterlab/kb_python@count-kite
```

    Collecting git+https://github.com/pachterlab/kb_python@count-kite
      Cloning https://github.com/pachterlab/kb_python (to revision count-kite) to /tmp/pip-req-build-ijboqjx7
      Running command git clone -q https://github.com/pachterlab/kb_python /tmp/pip-req-build-ijboqjx7
      Running command git checkout -b count-kite --track origin/count-kite
      Switched to a new branch 'count-kite'
      Branch 'count-kite' set up to track remote branch 'count-kite' from 'origin'.
    Collecting anndata>=0.6.22.post1
    [?25l  Downloading https://files.pythonhosted.org/packages/2b/72/87196c15f68d9865c31a43a10cf7c50bcbcedd5607d09f9aada0b3963103/anndata-0.6.22.post1-py3-none-any.whl (47kB)
    [K     |████████████████████████████████| 51kB 2.5MB/s 
    [?25hCollecting loompy>=3.0.6
    [?25l  Downloading https://files.pythonhosted.org/packages/36/52/74ed37ae5988522fbf87b856c67c4f80700e6452410b4cd80498c5f416f9/loompy-3.0.6.tar.gz (41kB)
    [K     |████████████████████████████████| 51kB 7.5MB/s 
    [?25hRequirement already satisfied: requests>=2.19.0 in /usr/local/lib/python3.6/dist-packages (from kb-python==0.24.4) (2.21.0)
    Collecting tqdm>=4.39.0
    [?25l  Downloading https://files.pythonhosted.org/packages/72/c9/7fc20feac72e79032a7c8138fd0d395dc6d8812b5b9edf53c3afd0b31017/tqdm-4.41.1-py2.py3-none-any.whl (56kB)
    [K     |████████████████████████████████| 61kB 8.0MB/s 
    [?25hRequirement already satisfied: natsort in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (5.5.0)
    Requirement already satisfied: numpy~=1.14 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (1.17.5)
    Requirement already satisfied: scipy~=1.0 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (1.4.1)
    Requirement already satisfied: pandas>=0.23.0 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (0.25.3)
    Requirement already satisfied: h5py in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (2.8.0)
    Requirement already satisfied: setuptools in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python==0.24.4) (42.0.2)
    Requirement already satisfied: numba in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python==0.24.4) (0.47.0)
    Requirement already satisfied: click in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python==0.24.4) (7.0)
    Collecting numpy-groupies
    [?25l  Downloading https://files.pythonhosted.org/packages/57/ae/18217b57ba3e4bb8a44ecbfc161ed065f6d1b90c75d404bd6ba8d6f024e2/numpy_groupies-0.9.10.tar.gz (43kB)
    [K     |████████████████████████████████| 51kB 6.7MB/s 
    [?25hRequirement already satisfied: chardet<3.1.0,>=3.0.2 in /usr/local/lib/python3.6/dist-packages (from requests>=2.19.0->kb-python==0.24.4) (3.0.4)
    Requirement already satisfied: certifi>=2017.4.17 in /usr/local/lib/python3.6/dist-packages (from requests>=2.19.0->kb-python==0.24.4) (2019.11.28)
    Requirement already satisfied: idna<2.9,>=2.5 in /usr/local/lib/python3.6/dist-packages (from requests>=2.19.0->kb-python==0.24.4) (2.8)
    Requirement already satisfied: urllib3<1.25,>=1.21.1 in /usr/local/lib/python3.6/dist-packages (from requests>=2.19.0->kb-python==0.24.4) (1.24.3)
    Requirement already satisfied: pytz>=2017.2 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python==0.24.4) (2018.9)
    Requirement already satisfied: python-dateutil>=2.6.1 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python==0.24.4) (2.6.1)
    Requirement already satisfied: six in /usr/local/lib/python3.6/dist-packages (from h5py->anndata>=0.6.22.post1->kb-python==0.24.4) (1.12.0)
    Requirement already satisfied: llvmlite>=0.31.0dev0 in /usr/local/lib/python3.6/dist-packages (from numba->loompy>=3.0.6->kb-python==0.24.4) (0.31.0)
    Building wheels for collected packages: kb-python, loompy, numpy-groupies
      Building wheel for kb-python (setup.py) ... [?25l[?25hdone
      Created wheel for kb-python: filename=kb_python-0.24.4-cp36-none-any.whl size=80381943 sha256=a321915d8605576908ac7e00a38e02a921d30f572ac59ed8208a7a69513ec758
      Stored in directory: /tmp/pip-ephem-wheel-cache-se6n68ix/wheels/8e/56/56/c89223de74af26792675e82f4bb5223e7cf0d653a33038e34c
      Building wheel for loompy (setup.py) ... [?25l[?25hdone
      Created wheel for loompy: filename=loompy-3.0.6-cp36-none-any.whl size=47896 sha256=f0f0f873d0617d273175b2894f232e39d81f169cbe525e86f502d3531aebf644
      Stored in directory: /root/.cache/pip/wheels/f9/a4/90/5a98ad83419732b0fba533b81a2a52ba3dbe230a936ca4cdc9
      Building wheel for numpy-groupies (setup.py) ... [?25l[?25hdone
      Created wheel for numpy-groupies: filename=numpy_groupies-0+unknown-cp36-none-any.whl size=28044 sha256=bfba2c0157a987e61dedbc3e60a7dee51422ae93e6b0bf173b0179514d761854
      Stored in directory: /root/.cache/pip/wheels/30/ac/83/64d5f9293aeaec63f9539142fc629a41af064cae1b3d8d94aa
    Successfully built kb-python loompy numpy-groupies
    Installing collected packages: anndata, numpy-groupies, loompy, tqdm, kb-python
      Found existing installation: tqdm 4.28.1
        Uninstalling tqdm-4.28.1:
          Successfully uninstalled tqdm-4.28.1
    Successfully installed anndata-0.6.22.post1 kb-python-0.24.4 loompy-3.0.6 numpy-groupies-0+unknown tqdm-4.41.1




### Build the feature barcode mismatch index

`kb` is able to generate a FASTA file containing all hamming distance < 2 variants of the feature barcodes and create a kallisto index of these sequences. But it in order to do so, we first need to prepare a TSV containing feature barcode sequences in the first column and the feature barcode names in the second.

First, we download the feature reference file provided by 10x Genomics.


```
!wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_feature_ref.csv
```

    --2020-01-15 23:32:15--  http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_feature_ref.csv
    Resolving cf.10xgenomics.com (cf.10xgenomics.com)... 13.224.29.124, 13.224.29.102, 13.224.29.56, ...
    Connecting to cf.10xgenomics.com (cf.10xgenomics.com)|13.224.29.124|:80... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 1473 (1.4K) [text/csv]
    Saving to: ‘pbmc_1k_protein_v3_feature_ref.csv’
    
    pbmc_1k_protein_v3_ 100%[===================>]   1.44K  --.-KB/s    in 0s      
    
    2020-01-15 23:32:15 (165 MB/s) - ‘pbmc_1k_protein_v3_feature_ref.csv’ saved [1473/1473]
    


Let's load it in as a Pandas DataFrame.


```
import pandas as pd

df = pd.read_csv('pbmc_1k_protein_v3_feature_ref.csv')
df
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
      <th>id</th>
      <th>name</th>
      <th>read</th>
      <th>pattern</th>
      <th>sequence</th>
      <th>feature_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>CD3</td>
      <td>CD3_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>AACAAGACCCTTGAG</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>1</th>
      <td>CD4</td>
      <td>CD4_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>TACCCGTAATAGCGT</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>2</th>
      <td>CD8a</td>
      <td>CD8a_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>ATTGGCACTCAGATG</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>3</th>
      <td>CD14</td>
      <td>CD14_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>GAAAGTCAAAGCACT</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>4</th>
      <td>CD15</td>
      <td>CD15_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>ACGAATCAATCTGTG</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>5</th>
      <td>CD16</td>
      <td>CD16_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>GTCTTTGTCAGTGCA</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>6</th>
      <td>CD56</td>
      <td>CD56_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>GTTGTCCGACAATAC</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>7</th>
      <td>CD19</td>
      <td>CD19_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>TCAACGCTTGGCTAG</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>8</th>
      <td>CD25</td>
      <td>CD25_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>GTGCATTCAACAGTA</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>9</th>
      <td>CD45RA</td>
      <td>CD45RA_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>GATGAGAACAGGTTT</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>10</th>
      <td>CD45RO</td>
      <td>CD45RO_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>TGCATGTCATCGGTG</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>11</th>
      <td>PD-1</td>
      <td>PD-1_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>AAGTCGTGAGGCATG</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>12</th>
      <td>TIGIT</td>
      <td>TIGIT_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>TGAAGGCTCATTTGT</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>13</th>
      <td>CD127</td>
      <td>CD127_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>ACATTGACGCAACTA</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>14</th>
      <td>IgG2a</td>
      <td>IgG2a_control_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>CTCTATTCAGACCAG</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>15</th>
      <td>IgG1</td>
      <td>IgG1_control_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>ACTCACTGGAGTCTC</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>16</th>
      <td>IgG2b</td>
      <td>IgG2b_control_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>ATCACATCGTTGCCA</td>
      <td>Antibody Capture</td>
    </tr>
  </tbody>
</table>
</div>



We'll convert this dataframe into a TSV format that `kb` requires.


```
df[['sequence', 'id']].to_csv('features.tsv', index=None, header=None, sep='\t')
!cat features.tsv
```

    AACAAGACCCTTGAG	CD3
    TACCCGTAATAGCGT	CD4
    ATTGGCACTCAGATG	CD8a
    GAAAGTCAAAGCACT	CD14
    ACGAATCAATCTGTG	CD15
    GTCTTTGTCAGTGCA	CD16
    GTTGTCCGACAATAC	CD56
    TCAACGCTTGGCTAG	CD19
    GTGCATTCAACAGTA	CD25
    GATGAGAACAGGTTT	CD45RA
    TGCATGTCATCGGTG	CD45RO
    AAGTCGTGAGGCATG	PD-1
    TGAAGGCTCATTTGT	TIGIT
    ACATTGACGCAACTA	CD127
    CTCTATTCAGACCAG	IgG2a
    ACTCACTGGAGTCTC	IgG1
    ATCACATCGTTGCCA	IgG2b


Finally, we use `kb` to generate the mismatch kallisto index.


```
!kb ref -i mismatch.idx -f1 mismatch.fa -g t2g.txt --workflow kite features.tsv
```

    [2020-01-15 23:32:21,669]    INFO Generating mismatch FASTA at mismatch.fa
    [2020-01-15 23:32:21,680]    INFO Creating transcript-to-gene mapping at t2g.txt
    [2020-01-15 23:32:21,684]    INFO Indexing to mismatch.idx


### Generate a feature count matrix in H5AD format

The following command will generate an RNA count matrix of cells (rows) by genes (columns) in H5AD format, which is a binary format used to store [Anndata](https://anndata.readthedocs.io/en/stable/) objects. Notice we are providing the index and transcript-to-gene mapping we generated in the previous step to the `-i` and `-g` arguments respectively. Also, these reads were generated with the 10x Genomics Chromium Single Cell v3 Chemistry, hence the `-x 10xv3` argument. To view other supported technologies, run `kb --list`.

__Note:__ If you would like a Loom file instead, replace the `--h5ad` flag with `--loom`. If you want to use the raw matrix output by `kb` instead of their H5AD or Loom converted files, omit these flags.


```
%%time
!kb count --h5ad -i mismatch.idx -g t2g.txt -x 10xv3 --workflow kite -t 2 \
pbmc_1k_protein_v3_antibody_S2_L001_R1_001.fastq.gz \
pbmc_1k_protein_v3_antibody_S2_L001_R2_001.fastq.gz \
pbmc_1k_protein_v3_antibody_S2_L002_R1_001.fastq.gz \
pbmc_1k_protein_v3_antibody_S2_L002_R2_001.fastq.gz
```

    [2020-01-15 23:32:28,148]    INFO Generating BUS file from
    [2020-01-15 23:32:28,148]    INFO         pbmc_1k_protein_v3_antibody_S2_L001_R1_001.fastq.gz
    [2020-01-15 23:32:28,148]    INFO         pbmc_1k_protein_v3_antibody_S2_L001_R2_001.fastq.gz
    [2020-01-15 23:32:28,148]    INFO         pbmc_1k_protein_v3_antibody_S2_L002_R1_001.fastq.gz
    [2020-01-15 23:32:28,148]    INFO         pbmc_1k_protein_v3_antibody_S2_L002_R2_001.fastq.gz
    [2020-01-15 23:34:23,527]    INFO Sorting BUS file ./output.bus to ./tmp/output.s.bus
    [2020-01-15 23:34:34,619]    INFO Whitelist not provided
    [2020-01-15 23:34:34,619]    INFO Copying pre-packaged 10XV3 whitelist to .
    [2020-01-15 23:34:35,535]    INFO Inspecting BUS file ./tmp/output.s.bus
    [2020-01-15 23:34:42,910]    INFO Correcting BUS records in ./tmp/output.s.bus to ./tmp/output.s.c.bus with whitelist ./10xv3_whitelist.txt
    [2020-01-15 23:34:59,040]    INFO Sorting BUS file ./tmp/output.s.c.bus to ./output.unfiltered.bus
    [2020-01-15 23:35:07,894]    INFO Generating count matrix ./counts_unfiltered/cells_x_features from BUS file ./output.unfiltered.bus
    [2020-01-15 23:35:10,231]    INFO Reading matrix ./counts_unfiltered/cells_x_features.mtx
    [2020-01-15 23:35:11,754]    INFO Writing matrix to h5ad ./counts_unfiltered/adata.h5ad
    CPU times: user 791 ms, sys: 93.8 ms, total: 884 ms
    Wall time: 2min 45s


## Analysis

In this part of the tutorial, we will load the RNA count matrix generated by `kb count` into Python and cluster the cells with Leiden.

### Install packages

Google Colab does not come with `Scanpy`, `python-igraph`, or `louvain` (but comes with `matplotlib`, `numpy`, `pandas`, and `scipy`).


```
!pip install leidenalg scanpy MulticoreTSNE
```

    Collecting leidenalg
    [?25l  Downloading https://files.pythonhosted.org/packages/b6/cc/d76baf78a3924ba6093a3ce8d14e2289f1d718bd3bcbb8252bb131d12daa/leidenalg-0.7.0.tar.gz (92kB)
    [K     |████████████████████████████████| 102kB 4.7MB/s 
    [?25hCollecting scanpy
    [?25l  Downloading https://files.pythonhosted.org/packages/ef/76/c3af48ffa7b657412d9b0cc8c13fa7918c41c26e89c45daf5aa8f3a475ee/scanpy-1.4.5.post2-py3-none-any.whl (310kB)
    [K     |████████████████████████████████| 317kB 15.5MB/s 
    [?25hCollecting MulticoreTSNE
      Downloading https://files.pythonhosted.org/packages/2d/e8/2afa896fa4eebfa1d0d0ba2673fddac45582ec0f06b2bdda88108ced5425/MulticoreTSNE-0.1.tar.gz
    Collecting python-igraph>=0.7.1.0
    [?25l  Downloading https://files.pythonhosted.org/packages/0f/a0/4e7134f803737aa6eebb4e5250565ace0e2599659e22be7f7eba520ff017/python-igraph-0.7.1.post6.tar.gz (377kB)
    [K     |████████████████████████████████| 378kB 41.3MB/s 
    [?25hRequirement already satisfied: pandas>=0.21 in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.25.3)
    Requirement already satisfied: tqdm in /usr/local/lib/python3.6/dist-packages (from scanpy) (4.41.1)
    Requirement already satisfied: scikit-learn>=0.21.2 in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.22.1)
    Requirement already satisfied: statsmodels>=0.10.0rc2 in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.10.2)
    Requirement already satisfied: patsy in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.5.1)
    Requirement already satisfied: numba>=0.41.0 in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.47.0)
    Requirement already satisfied: importlib-metadata>=0.7; python_version < "3.8" in /usr/local/lib/python3.6/dist-packages (from scanpy) (1.4.0)
    Collecting setuptools-scm
      Downloading https://files.pythonhosted.org/packages/1d/70/97966deebaeeda0b81d3cd63ba9f8ec929b838871ed17476de9d8159db3e/setuptools_scm-3.3.3-py2.py3-none-any.whl
    Requirement already satisfied: joblib in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.14.1)
    Collecting legacy-api-wrap
      Downloading https://files.pythonhosted.org/packages/a4/68/da997bc56bb69dcdcee4054f0bc42266909307b905389fbc54c9158f42da/legacy_api_wrap-1.2-py3-none-any.whl
    Collecting matplotlib==3.0.*
    [?25l  Downloading https://files.pythonhosted.org/packages/e9/69/f5e05f578585ed9935247be3788b374f90701296a70c8871bcd6d21edb00/matplotlib-3.0.3-cp36-cp36m-manylinux1_x86_64.whl (13.0MB)
    [K     |████████████████████████████████| 13.0MB 54.1MB/s 
    [?25hRequirement already satisfied: h5py!=2.10.0 in /usr/local/lib/python3.6/dist-packages (from scanpy) (2.8.0)
    Requirement already satisfied: seaborn in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.9.0)
    Requirement already satisfied: tables in /usr/local/lib/python3.6/dist-packages (from scanpy) (3.4.4)
    Requirement already satisfied: networkx in /usr/local/lib/python3.6/dist-packages (from scanpy) (2.4)
    Requirement already satisfied: scipy>=1.3 in /usr/local/lib/python3.6/dist-packages (from scanpy) (1.4.1)
    Requirement already satisfied: packaging in /usr/local/lib/python3.6/dist-packages (from scanpy) (20.0)
    Requirement already satisfied: umap-learn>=0.3.10 in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.3.10)
    Requirement already satisfied: anndata>=0.6.22.post1 in /usr/local/lib/python3.6/dist-packages (from scanpy) (0.6.22.post1)
    Requirement already satisfied: natsort in /usr/local/lib/python3.6/dist-packages (from scanpy) (5.5.0)
    Requirement already satisfied: numpy in /usr/local/lib/python3.6/dist-packages (from MulticoreTSNE) (1.17.5)
    Requirement already satisfied: cffi in /usr/local/lib/python3.6/dist-packages (from MulticoreTSNE) (1.13.2)
    Requirement already satisfied: python-dateutil>=2.6.1 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.21->scanpy) (2.6.1)
    Requirement already satisfied: pytz>=2017.2 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.21->scanpy) (2018.9)
    Requirement already satisfied: six in /usr/local/lib/python3.6/dist-packages (from patsy->scanpy) (1.12.0)
    Requirement already satisfied: llvmlite>=0.31.0dev0 in /usr/local/lib/python3.6/dist-packages (from numba>=0.41.0->scanpy) (0.31.0)
    Requirement already satisfied: setuptools in /usr/local/lib/python3.6/dist-packages (from numba>=0.41.0->scanpy) (42.0.2)
    Requirement already satisfied: zipp>=0.5 in /usr/local/lib/python3.6/dist-packages (from importlib-metadata>=0.7; python_version < "3.8"->scanpy) (0.6.0)
    Collecting get-version>=2.0.4
    [?25l  Downloading https://files.pythonhosted.org/packages/23/48/7610e884e62fff2183e7bc8592397c39a020267fb5147905fcd3f9cc820c/get_version-2.1-py3-none-any.whl (43kB)
    [K     |████████████████████████████████| 51kB 7.0MB/s 
    [?25hRequirement already satisfied: cycler>=0.10 in /usr/local/lib/python3.6/dist-packages (from matplotlib==3.0.*->scanpy) (0.10.0)
    Requirement already satisfied: pyparsing!=2.0.4,!=2.1.2,!=2.1.6,>=2.0.1 in /usr/local/lib/python3.6/dist-packages (from matplotlib==3.0.*->scanpy) (2.4.6)
    Requirement already satisfied: kiwisolver>=1.0.1 in /usr/local/lib/python3.6/dist-packages (from matplotlib==3.0.*->scanpy) (1.1.0)
    Requirement already satisfied: numexpr>=2.5.2 in /usr/local/lib/python3.6/dist-packages (from tables->scanpy) (2.7.1)
    Requirement already satisfied: decorator>=4.3.0 in /usr/local/lib/python3.6/dist-packages (from networkx->scanpy) (4.4.1)
    Requirement already satisfied: pycparser in /usr/local/lib/python3.6/dist-packages (from cffi->MulticoreTSNE) (2.19)
    Requirement already satisfied: more-itertools in /usr/local/lib/python3.6/dist-packages (from zipp>=0.5->importlib-metadata>=0.7; python_version < "3.8"->scanpy) (8.0.2)
    Building wheels for collected packages: leidenalg, MulticoreTSNE, python-igraph
      Building wheel for leidenalg (setup.py) ... [?25l[?25hdone
      Created wheel for leidenalg: filename=leidenalg-0.7.0-cp36-cp36m-linux_x86_64.whl size=1104867 sha256=9fa39f77fd8a43de76be0788ec28b1273681be13da29dfe558c876eead3be36f
      Stored in directory: /root/.cache/pip/wheels/29/55/48/5a04693a10f50297bcda23819ca23ab3470a61dd911851c8bd
      Building wheel for MulticoreTSNE (setup.py) ... [?25l[?25hdone
      Created wheel for MulticoreTSNE: filename=MulticoreTSNE-0.1-cp36-cp36m-linux_x86_64.whl size=68507 sha256=05e5da74abe7e90b95993cc4d7cf19f320cc84225dd53b5efc94cbc733004b5c
      Stored in directory: /root/.cache/pip/wheels/27/59/53/3b52ee63add3692254c30d687fa4dff4d128d0557861fb028e
      Building wheel for python-igraph (setup.py) ... [?25l[?25hdone
      Created wheel for python-igraph: filename=python_igraph-0.7.1.post6-cp36-cp36m-linux_x86_64.whl size=2224304 sha256=cc023e1a4a81a0938a8bf55f988c4b4bacf24932611019f23de8def21ba8e0f7
      Stored in directory: /root/.cache/pip/wheels/41/d6/02/34eebae97e25f5b87d60f4c0687e00523e3f244fa41bc3f4a7
    Successfully built leidenalg MulticoreTSNE python-igraph
    [31mERROR: plotnine 0.6.0 has requirement matplotlib>=3.1.1, but you'll have matplotlib 3.0.3 which is incompatible.[0m
    [31mERROR: mizani 0.6.0 has requirement matplotlib>=3.1.1, but you'll have matplotlib 3.0.3 which is incompatible.[0m
    [31mERROR: albumentations 0.1.12 has requirement imgaug<0.2.7,>=0.2.5, but you'll have imgaug 0.2.9 which is incompatible.[0m
    Installing collected packages: python-igraph, leidenalg, setuptools-scm, get-version, legacy-api-wrap, matplotlib, scanpy, MulticoreTSNE
      Found existing installation: matplotlib 3.1.2
        Uninstalling matplotlib-3.1.2:
          Successfully uninstalled matplotlib-3.1.2
    Successfully installed MulticoreTSNE-0.1 get-version-2.1 legacy-api-wrap-1.2 leidenalg-0.7.0 matplotlib-3.0.3 python-igraph-0.7.1.post6 scanpy-1.4.5.post2 setuptools-scm-3.3.3




### Import packages


```
import anndata
import numpy as np
import scanpy as sc
```


```
adata = anndata.read_h5ad('counts_unfiltered/adata.h5ad')
```


```
adata
```




    AnnData object with n_obs × n_vars = 124716 × 17 




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
    </tr>
    <tr>
      <th>index</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>AAACCCAAGAAACCCA</th>
    </tr>
    <tr>
      <th>AAACCCAAGACGAGGA</th>
    </tr>
    <tr>
      <th>AAACCCAAGAGTGTGT</th>
    </tr>
    <tr>
      <th>AAACCCAAGAGTGTTG</th>
    </tr>
    <tr>
      <th>AAACCCAAGATAGCAC</th>
    </tr>
  </tbody>
</table>
</div>




```
adata.var
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
    <tr>
      <th>index</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>CD3</th>
    </tr>
    <tr>
      <th>CD4</th>
    </tr>
    <tr>
      <th>CD8a</th>
    </tr>
    <tr>
      <th>CD14</th>
    </tr>
    <tr>
      <th>CD15</th>
    </tr>
    <tr>
      <th>CD16</th>
    </tr>
    <tr>
      <th>CD56</th>
    </tr>
    <tr>
      <th>CD19</th>
    </tr>
    <tr>
      <th>CD25</th>
    </tr>
    <tr>
      <th>CD45RA</th>
    </tr>
    <tr>
      <th>CD45RO</th>
    </tr>
    <tr>
      <th>PD-1</th>
    </tr>
    <tr>
      <th>TIGIT</th>
    </tr>
    <tr>
      <th>CD127</th>
    </tr>
    <tr>
      <th>IgG2a</th>
    </tr>
    <tr>
      <th>IgG1</th>
    </tr>
    <tr>
      <th>IgG2b</th>
    </tr>
  </tbody>
</table>
</div>



### Plot counts


```
sc.pp.filter_cells(adata, min_counts=0)
```


```
sc.pp.filter_genes(adata, min_counts=0)
```


```
sc.pl.violin(adata, keys='n_counts')
```


![png](kb_kite_files/kb_kite_31_0.png)



```
adata.obs['n_countslog'] = np.log1p(adata.obs['n_counts'])
```


```
sc.pl.violin(adata, keys='n_countslog')
```


![png](kb_kite_files/kb_kite_33_0.png)



```
adata.obs.index
```




    Index(['AAACCCAAGAAACCCA', 'AAACCCAAGACGAGGA', 'AAACCCAAGAGTGTGT',
           'AAACCCAAGAGTGTTG', 'AAACCCAAGATAGCAC', 'AAACCCAAGATGAGTC',
           'AAACCCAAGATGGTAC', 'AAACCCAAGATTCGTT', 'AAACCCAAGATTTGGG',
           'AAACCCAAGCAAGCAT',
           ...
           'TTTGTTGTCGTCGCCT', 'TTTGTTGTCGTTGACG', 'TTTGTTGTCTAACCGG',
           'TTTGTTGTCTATGTAG', 'TTTGTTGTCTCAACAA', 'TTTGTTGTCTCACTCA',
           'TTTGTTGTCTCTTCGA', 'TTTGTTGTCTCTTGGT', 'TTTGTTGTCTGCACTT',
           'TTTGTTGTCTGCGACA'],
          dtype='object', name='index', length=124716)




```
sc.pp.filter_cells(adata, min_counts=1000)
sc.pl.violin(adata, keys='n_countslog', title="kallisto UMI counts")
adata
```


![png](kb_kite_files/kb_kite_35_0.png)





    AnnData object with n_obs × n_vars = 725 × 17 
        obs: 'n_counts', 'n_countslog'
        var: 'n_counts'



Here are violin plots for each Feature Barcode (antibody-oligo conjugates, x-axis) across all cells.


```
sc.pl.violin(adata, keys=list(adata.var.index)[-17:], xlabel='kallisto')
```


![png](kb_kite_files/kb_kite_37_0.png)


### Cluster with Leiden


```
sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000)
```


```
sc.pp.neighbors(adata)
sc.tl.umap(adata)
```


```
sc.tl.leiden(adata, resolution=0.05)
```


```
sc.pl.umap(adata, color='leiden', palette='tab10')
```


![png](kb_kite_files/kb_kite_42_0.png)


### Embedding and Antibody Quantification


```
sc.pl.umap(adata, color=adata.var.index)
```


![png](kb_kite_files/kb_kite_44_0.png)



```
sc.pl.violin(adata, keys=list(adata.var.index[:2]), groupby='leiden')
```


![png](kb_kite_files/kb_kite_45_0.png)



```

```
