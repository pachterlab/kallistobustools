# Pre-processing mouse single-nuclei RNA-seq data with kallisto and bustools

In this tutorial we will process the 10x dataset [1k Brain Nuclei from an E18 Mouse](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/nuclei_900) using kallisto bus and a custom built DNA and intron index for mouse. We will generate two matrices: one for spliced transcripts and one for unspliced transcripts, and sum them to obtain total nuclear transcripts.


```
!date
```

    Thu Jan 16 15:52:05 UTC 2020


## Pre-processing

### Download the data

__Note:__ We use the `-O` option for `wget` to rename the files to easily identify them.


```
%%time
!wget https://caltech.box.com/shared/static/j337aflq9ublmwaripkepob41mr23216.txt -O checksums.txt
!wget https://caltech.box.com/shared/static/2j8shgwmalzcjawuow51678a8yssvdef.gz -O nuclei_900_S1_L001_R1_001.fastq.gz
!wget https://caltech.box.com/shared/static/k2yydqlz2jtckw1shk5h536mxn47nm9n.gz -O nuclei_900_S1_L001_R2_001.fastq.gz
!wget https://caltech.box.com/shared/static/tlqdm0w3tvy8ogyktsz7ahggwurc6kkj.gz -O nuclei_900_S1_L002_R1_001.fastq.gz
!wget https://caltech.box.com/shared/static/gqrvkqllr9d7zq4e3yfrng9kgfbejowe.gz -O nuclei_900_S1_L002_R2_001.fastq.gz
```

    --2020-01-16 15:52:07--  https://caltech.box.com/shared/static/j337aflq9ublmwaripkepob41mr23216.txt
    Resolving caltech.box.com (caltech.box.com)... 107.152.26.197, 107.152.27.197
    Connecting to caltech.box.com (caltech.box.com)|107.152.26.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/j337aflq9ublmwaripkepob41mr23216.txt [following]
    --2020-01-16 15:52:12--  https://caltech.box.com/public/static/j337aflq9ublmwaripkepob41mr23216.txt
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/j337aflq9ublmwaripkepob41mr23216.txt [following]
    --2020-01-16 15:52:12--  https://caltech.app.box.com/public/static/j337aflq9ublmwaripkepob41mr23216.txt
    Resolving caltech.app.box.com (caltech.app.box.com)... 107.152.27.199, 107.152.26.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|107.152.27.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!GcwiYMLuutOlK9x2GOcgBeJEzw7h5lJWU6P_DhY5cEFRktRc-MWgWWB1g94XLVSAuhxsYnAcFDQJ4VBAH_a3ObNnzvrJJwFkt220Nma3TJoLifXB3oWGWH3teBvR3pcjKZ7KijVqs6vO7Yofh44COuvzpbvEYfUqaXK-8PqHHQtQ7ZavLSV20ykoB81K-ZyLzSOIFoJjG8o-D8uz3uw5AInnVlDuKdN1qyQ17eujvTTBn4sj4joxaifmiZK-O7XA-FbcXkdHYUKwNhWzLPexMBsJKd1WT5qF7VAVmhajA6wo46bsrF_LTGETw_U5difW6etQNE_IuthI4anrdPV4ClJQDGCrL8whiuR5EGttUDzBmdNXgKJ8fTgOMUVQ5y4Ceij-F4bXTUW9lGGSYok7gdyHt3VnjkdMCbvq-t9kPN8esW_deFaJAH3BzU_jZ6u4652y5MDrWZMIg6r7wxMpGkG0smbZHN0PDIkESjMmkgSC-3mDxucly0sxtcD5-er0lOVpPuJivePQl_rNBSKmbXOrpSNfkWcR-o6eFkXVQcdHw-1j_DUpBtOqfhIg0rVCkvI82UBKesxf0gbVBERcQ54RTRdfDqfuXZKwbogxTie85n5LC3F75QVtJr9TT8y18nEgj02qEoChVIcIBJiiTv7pTAW_rCxBlOGDNN-iY4VyeP2PBZ6CKkEHlbaFT2L0LYLQMfvkXBYEC00S0e2rqQoXED6sTvxVVD3Df4tMskd9l7NDIF1OxNk7IueJhtiDICKMydwA860XmmOkfGMK1YkYJKdpqoAU04nY2_SO0gmtjEtJMSqxcr6Yr2nuVbnmU1sRu0Su2lMIrIB1QLF7ma1gsYArHbBa3e8tGddxykgkA_8gfUbRFzmmqvKKCWYo1yPnGqoZUNmg-s26zQimtE-ahQ4ePIA35N9HCvT3SP2IPIr1lLa1zkKckGFLSRvrlCV0X5Hr-I2CJlgslvnER35EIIxnnKA-Ey7wSrMbAOBkBVBMLpLGL15on0e-tP1knNZ2JdqQ0WahXYJaAmOk5ks4Dar8pO0gwsQIibD5BmtTW6Ag9qg0Au0Zgg5mqWUJc-QTsh9eqAx73I-YXBAci_oaKkmv4uZjej2NHM0lq8mej4onFQNAsA_m-Zxvgx5tdv-FqfcVvPSOMw10_wWatDcXxF-TcjTUtDpqhIJDfJNq8_CPzX_XyW6j05F5-Qa-pCZZjAHjPL8K0Mtnvj8Y2_2NIJW7kxyerJ2YUM_UipTYzxXVEV8wv4bCF5OxwV8fNlXnYkbl-2hc2NjTD9of0XA6uWa_FyPUetjM1xtQN89D3olZ9wyX411FkvCnj2qhBEs1vJjjJs4dT3pEY0MzO0CUtAgxOKppttE./download [following]
    --2020-01-16 15:52:13--  https://public.boxcloud.com/d/1/b1!GcwiYMLuutOlK9x2GOcgBeJEzw7h5lJWU6P_DhY5cEFRktRc-MWgWWB1g94XLVSAuhxsYnAcFDQJ4VBAH_a3ObNnzvrJJwFkt220Nma3TJoLifXB3oWGWH3teBvR3pcjKZ7KijVqs6vO7Yofh44COuvzpbvEYfUqaXK-8PqHHQtQ7ZavLSV20ykoB81K-ZyLzSOIFoJjG8o-D8uz3uw5AInnVlDuKdN1qyQ17eujvTTBn4sj4joxaifmiZK-O7XA-FbcXkdHYUKwNhWzLPexMBsJKd1WT5qF7VAVmhajA6wo46bsrF_LTGETw_U5difW6etQNE_IuthI4anrdPV4ClJQDGCrL8whiuR5EGttUDzBmdNXgKJ8fTgOMUVQ5y4Ceij-F4bXTUW9lGGSYok7gdyHt3VnjkdMCbvq-t9kPN8esW_deFaJAH3BzU_jZ6u4652y5MDrWZMIg6r7wxMpGkG0smbZHN0PDIkESjMmkgSC-3mDxucly0sxtcD5-er0lOVpPuJivePQl_rNBSKmbXOrpSNfkWcR-o6eFkXVQcdHw-1j_DUpBtOqfhIg0rVCkvI82UBKesxf0gbVBERcQ54RTRdfDqfuXZKwbogxTie85n5LC3F75QVtJr9TT8y18nEgj02qEoChVIcIBJiiTv7pTAW_rCxBlOGDNN-iY4VyeP2PBZ6CKkEHlbaFT2L0LYLQMfvkXBYEC00S0e2rqQoXED6sTvxVVD3Df4tMskd9l7NDIF1OxNk7IueJhtiDICKMydwA860XmmOkfGMK1YkYJKdpqoAU04nY2_SO0gmtjEtJMSqxcr6Yr2nuVbnmU1sRu0Su2lMIrIB1QLF7ma1gsYArHbBa3e8tGddxykgkA_8gfUbRFzmmqvKKCWYo1yPnGqoZUNmg-s26zQimtE-ahQ4ePIA35N9HCvT3SP2IPIr1lLa1zkKckGFLSRvrlCV0X5Hr-I2CJlgslvnER35EIIxnnKA-Ey7wSrMbAOBkBVBMLpLGL15on0e-tP1knNZ2JdqQ0WahXYJaAmOk5ks4Dar8pO0gwsQIibD5BmtTW6Ag9qg0Au0Zgg5mqWUJc-QTsh9eqAx73I-YXBAci_oaKkmv4uZjej2NHM0lq8mej4onFQNAsA_m-Zxvgx5tdv-FqfcVvPSOMw10_wWatDcXxF-TcjTUtDpqhIJDfJNq8_CPzX_XyW6j05F5-Qa-pCZZjAHjPL8K0Mtnvj8Y2_2NIJW7kxyerJ2YUM_UipTYzxXVEV8wv4bCF5OxwV8fNlXnYkbl-2hc2NjTD9of0XA6uWa_FyPUetjM1xtQN89D3olZ9wyX411FkvCnj2qhBEs1vJjjJs4dT3pEY0MzO0CUtAgxOKppttE./download
    Resolving public.boxcloud.com (public.boxcloud.com)... 107.152.26.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|107.152.26.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 408 [text/plain]
    Saving to: ‘checksums.txt’
    
    checksums.txt       100%[===================>]     408  --.-KB/s    in 0s      
    
    2020-01-16 15:52:13 (95.3 MB/s) - ‘checksums.txt’ saved [408/408]
    
    --2020-01-16 15:52:14--  https://caltech.box.com/shared/static/2j8shgwmalzcjawuow51678a8yssvdef.gz
    Resolving caltech.box.com (caltech.box.com)... 107.152.26.197, 107.152.27.197
    Connecting to caltech.box.com (caltech.box.com)|107.152.26.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/2j8shgwmalzcjawuow51678a8yssvdef.gz [following]
    --2020-01-16 15:52:14--  https://caltech.box.com/public/static/2j8shgwmalzcjawuow51678a8yssvdef.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/2j8shgwmalzcjawuow51678a8yssvdef.gz [following]
    --2020-01-16 15:52:14--  https://caltech.app.box.com/public/static/2j8shgwmalzcjawuow51678a8yssvdef.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 107.152.27.199, 107.152.26.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|107.152.27.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!5U6Zu_0JeqnVVcz5PDt4XxrXXLpaSp2XNMDZ8Xf3Z1OJcH_gaUEBqmCnFYkZmilbCIlf_R01pEelr0aDWUVm5FliiDmW4jTzQnM_hm2iEesumpJG8ssdu40dAtXERs9ZHNbHSy1P9nXJGtccGJPrKBCQBQyBBCjHigChYvrcrxvArQBzjatym5JkG2cZagUV6Dn8VA4N9LCShceA3YPhlXpG_-X7p1T1DcIoi6U6jbZUMbJW_NQWPvyELi1sfgy4_WL3ubWZnTDroUlDAAGlJ2A38LW9pljkIN2tdupF3yuvulv2MPfz0Bn_K5ncwoicKNNZ95k3CELi_70hkWmcHVvtJEO0CvI6r8rHFdghxsUFxwopNq-77sNYnxEZiCE-lj7bKZzzhH015p7Kx6haSmbe9WqmVkJMP1Cn1hke_DTCMv_d13BJkvt5Ywt1uPf0n-k3bmOUZpoAa4fXTWLrE0fn4PtmgXRZ4RZURHcIT1WCO8FKRnkUhlp-Vg0W9n3pEwJUUivqO_tcav6KPYWVElVk67teYNNnEAeNF_8ogUvuITdrZWTWSpZ5agP-jFND4_1vf_L7ZvrOFwCgIz0SRByJdT9wCJ9ofGCyQ_ECuIpAiArn5BGvnNwutuuH3N3RgWOVTr5Mhy6zYbME5u1qocmUx0Uvv8CNnASSMQ74nMqLVoDd7EByDCrW9ONOn6MEW-tX4U65ArB3tnwT_Gnx7yvOgWNpIg8bjlb7FPX3kbvO7wQe7qMQqH0GvFveeQw3-3g04tya_WLMpcrd2dxdSfku3vkWeOd34jIL2IKVjd7wqUhf6O2DFQcRFrrZ7cxgiRXQ7uABbY7qS9eykxQvetzIhHbGSZEm9d0VSEQKEIb-KMimKxGVCul2FQSf-6beI6n-DBwso3VO7sxUzM8OCUN0YElsYuqfXVT2-5U_Hr-qkYQWi1S6xhGdif5toSuTEzn4_WPGTklAkkDHJT_Z83GMhLs65QCOfrELlsm_Tw0LWw-CpbRrrbAGbrZRCxDdh4NuHmV5aH0mREQgxQ1gta-9hrlBTTELfX73t27xgOBeAe2AazG6K8G8pDxxchRWUOjNHEwYJQ-rMT9n6xRHHJzm1NiOElOXU-AiI2012vx02xbJGo8vm3SQfV4On09G7Sblq_4D9eRHLzHf0dGiURAMrpiApc1VFC0-NlEhEBO-3TkiM_xI9OjZYYx9uiE9tEPmLS9zJUrnvoIbf5iHnrEha4-CFZw-KIUmYdlpTlSSILnT5mT1wzm_DUoMoIe9GlrGCRvx60cm6F5F2bC2ehrUqtJoeOZP9XCp7wW3V8T5vOLBgiUECr1Ry9hlZaqkFwWYDNHxwUabYmHNyAfTDpUI4zvF79mZ7b7sb2h_UpSmS4Ab_HLroZYBCe7iLn3qRuLoP_A3BvohWxHAO_dMYqV1VA../download [following]
    --2020-01-16 15:52:15--  https://public.boxcloud.com/d/1/b1!5U6Zu_0JeqnVVcz5PDt4XxrXXLpaSp2XNMDZ8Xf3Z1OJcH_gaUEBqmCnFYkZmilbCIlf_R01pEelr0aDWUVm5FliiDmW4jTzQnM_hm2iEesumpJG8ssdu40dAtXERs9ZHNbHSy1P9nXJGtccGJPrKBCQBQyBBCjHigChYvrcrxvArQBzjatym5JkG2cZagUV6Dn8VA4N9LCShceA3YPhlXpG_-X7p1T1DcIoi6U6jbZUMbJW_NQWPvyELi1sfgy4_WL3ubWZnTDroUlDAAGlJ2A38LW9pljkIN2tdupF3yuvulv2MPfz0Bn_K5ncwoicKNNZ95k3CELi_70hkWmcHVvtJEO0CvI6r8rHFdghxsUFxwopNq-77sNYnxEZiCE-lj7bKZzzhH015p7Kx6haSmbe9WqmVkJMP1Cn1hke_DTCMv_d13BJkvt5Ywt1uPf0n-k3bmOUZpoAa4fXTWLrE0fn4PtmgXRZ4RZURHcIT1WCO8FKRnkUhlp-Vg0W9n3pEwJUUivqO_tcav6KPYWVElVk67teYNNnEAeNF_8ogUvuITdrZWTWSpZ5agP-jFND4_1vf_L7ZvrOFwCgIz0SRByJdT9wCJ9ofGCyQ_ECuIpAiArn5BGvnNwutuuH3N3RgWOVTr5Mhy6zYbME5u1qocmUx0Uvv8CNnASSMQ74nMqLVoDd7EByDCrW9ONOn6MEW-tX4U65ArB3tnwT_Gnx7yvOgWNpIg8bjlb7FPX3kbvO7wQe7qMQqH0GvFveeQw3-3g04tya_WLMpcrd2dxdSfku3vkWeOd34jIL2IKVjd7wqUhf6O2DFQcRFrrZ7cxgiRXQ7uABbY7qS9eykxQvetzIhHbGSZEm9d0VSEQKEIb-KMimKxGVCul2FQSf-6beI6n-DBwso3VO7sxUzM8OCUN0YElsYuqfXVT2-5U_Hr-qkYQWi1S6xhGdif5toSuTEzn4_WPGTklAkkDHJT_Z83GMhLs65QCOfrELlsm_Tw0LWw-CpbRrrbAGbrZRCxDdh4NuHmV5aH0mREQgxQ1gta-9hrlBTTELfX73t27xgOBeAe2AazG6K8G8pDxxchRWUOjNHEwYJQ-rMT9n6xRHHJzm1NiOElOXU-AiI2012vx02xbJGo8vm3SQfV4On09G7Sblq_4D9eRHLzHf0dGiURAMrpiApc1VFC0-NlEhEBO-3TkiM_xI9OjZYYx9uiE9tEPmLS9zJUrnvoIbf5iHnrEha4-CFZw-KIUmYdlpTlSSILnT5mT1wzm_DUoMoIe9GlrGCRvx60cm6F5F2bC2ehrUqtJoeOZP9XCp7wW3V8T5vOLBgiUECr1Ry9hlZaqkFwWYDNHxwUabYmHNyAfTDpUI4zvF79mZ7b7sb2h_UpSmS4Ab_HLroZYBCe7iLn3qRuLoP_A3BvohWxHAO_dMYqV1VA../download
    Resolving public.boxcloud.com (public.boxcloud.com)... 107.152.26.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|107.152.26.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 627431870 (598M) [application/octet-stream]
    Saving to: ‘nuclei_900_S1_L001_R1_001.fastq.gz’
    
    nuclei_900_S1_L001_ 100%[===================>] 598.37M  23.1MB/s    in 26s     
    
    2020-01-16 15:52:41 (22.7 MB/s) - ‘nuclei_900_S1_L001_R1_001.fastq.gz’ saved [627431870/627431870]
    
    --2020-01-16 15:52:42--  https://caltech.box.com/shared/static/k2yydqlz2jtckw1shk5h536mxn47nm9n.gz
    Resolving caltech.box.com (caltech.box.com)... 107.152.26.197, 107.152.27.197
    Connecting to caltech.box.com (caltech.box.com)|107.152.26.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/k2yydqlz2jtckw1shk5h536mxn47nm9n.gz [following]
    --2020-01-16 15:52:42--  https://caltech.box.com/public/static/k2yydqlz2jtckw1shk5h536mxn47nm9n.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/k2yydqlz2jtckw1shk5h536mxn47nm9n.gz [following]
    --2020-01-16 15:52:42--  https://caltech.app.box.com/public/static/k2yydqlz2jtckw1shk5h536mxn47nm9n.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 107.152.27.199, 107.152.26.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|107.152.27.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!cM7F-H9GC5xcOTIcPEDIng1keb1_sf_4S9a_OTqkd7FRTT5sHg-gF9tD7JUuRRkJrhgOgAInFK8YcJQq_XA619I9Qj0UFxrLb7nMBToxBXF4Ml_gXYKiuq0Y6L2B65mL6mN1hROJDOl6X598o2N4nurQrUV5Hd2oABh70YznMfdDgX-U_Y-dQBMWL9EzLL2QPiYskaeHcNk1QdxoiciFCd6PwpXct-aIFfL1tHnlE9pJ2Jniad7y1Ymv8rdAN_nvzGq68-vhYAR6Ypsy3UeJ_84tMCMy2KdvY0aqQdnl7jpPLSHNDV1VWFLCJat2mYmwi4mx9jOI__ydx_efWv_T5MegHs3e1GmCA38YvswdJDjgqiChEuXwV6yvVmRZSVisScRJz5_LwwMKTmGmMu4Hpa0WsG0Mret8kZ8lujx4EZLHE4LQhUQ1aBO_xwDyPrZ_TavMDuAF_e2ZDkTahwALIF2cDVOYy58wzpa-iy7rbA5E-jlceLKd5IoKbanw0OHNRBE8Gi1pQ9-uP-rfg8wf3DU7qHh-C7-I-IJFDpmlfqNESMHBywNZmuOvcQ0hT6lDGMU0NZD4Y8WFjH9UmKE2apmwv5qyNFFXpyYsG2ne-zdLoe0FirKHDkto7hnX-3jewG9eY5V2oIMby7XBmJa2CQFL-Fz1VCiuD0HkkUobMFKzPAIsNpZJlUwfpNq21dJG9WLRzcQ-njX6A_hD621Zg_QGxt3kDZF_r7NpAfwSLz4aHyVjuBucQgbEVwiMtm-qUC5dsGJEv8qioS9dSC3JHuoYi7yde28uPW3RUGgDzD2ZDspMaEI7GQqrW_kqNHq3uTTG_hdx_o8e63dFZ0HXjWT-D09ziZd7_6A2p4NCf2aDNW_Ys8y0PkrVX1VyLHxqpT2j62-t9y5l-fd3AlRkgAe1WnmJVOUTKATHzYKK4MQKSW5bQ6fMV1Khb_HjVWzk-zpoZGBKNYkPEAGmQQv7Tmr9CZ_gD2aoCDibXCFBXi2k-TftSsnlzKbZkIF2IdSPHWpthXwj6Xf44Qnx1HRCp10MlrZmPI49TlNa2ZzHQUm2TEvFuqwzUMB7qIPWUEGdgSfIUpc00t3Yhx_2BaGbIzCp6Cp4O6BgTymoRSErj_Xg29Y4Xzkyen1nUc2zc0pWSd_h-iK5Jt9FmOupM4LJVSM_b04nIPARZFR7GxhRbzE5YX-EZhy2_yJIsH_Ajl0fPdZ2sDm3zne3OhdcBZ-Ea0BEqpemjOh1a717iX1b2LbLf5wkx2lmlVUhkO09FtwYqw-U26jrElqmfOl6vEDvG-GUJoYTHrSJ_nVFbaqt2eGqU6HGSyS7FQa0LVR5jANJPU2vi78s8y0FF75_CIMYHpjVYyQXkgEGNVmpm5wKBdagk4E9YSeSRey0gETkz1iVigAf80EL-c5nvMUfCvOktgDrjRo./download [following]
    --2020-01-16 15:52:43--  https://public.boxcloud.com/d/1/b1!cM7F-H9GC5xcOTIcPEDIng1keb1_sf_4S9a_OTqkd7FRTT5sHg-gF9tD7JUuRRkJrhgOgAInFK8YcJQq_XA619I9Qj0UFxrLb7nMBToxBXF4Ml_gXYKiuq0Y6L2B65mL6mN1hROJDOl6X598o2N4nurQrUV5Hd2oABh70YznMfdDgX-U_Y-dQBMWL9EzLL2QPiYskaeHcNk1QdxoiciFCd6PwpXct-aIFfL1tHnlE9pJ2Jniad7y1Ymv8rdAN_nvzGq68-vhYAR6Ypsy3UeJ_84tMCMy2KdvY0aqQdnl7jpPLSHNDV1VWFLCJat2mYmwi4mx9jOI__ydx_efWv_T5MegHs3e1GmCA38YvswdJDjgqiChEuXwV6yvVmRZSVisScRJz5_LwwMKTmGmMu4Hpa0WsG0Mret8kZ8lujx4EZLHE4LQhUQ1aBO_xwDyPrZ_TavMDuAF_e2ZDkTahwALIF2cDVOYy58wzpa-iy7rbA5E-jlceLKd5IoKbanw0OHNRBE8Gi1pQ9-uP-rfg8wf3DU7qHh-C7-I-IJFDpmlfqNESMHBywNZmuOvcQ0hT6lDGMU0NZD4Y8WFjH9UmKE2apmwv5qyNFFXpyYsG2ne-zdLoe0FirKHDkto7hnX-3jewG9eY5V2oIMby7XBmJa2CQFL-Fz1VCiuD0HkkUobMFKzPAIsNpZJlUwfpNq21dJG9WLRzcQ-njX6A_hD621Zg_QGxt3kDZF_r7NpAfwSLz4aHyVjuBucQgbEVwiMtm-qUC5dsGJEv8qioS9dSC3JHuoYi7yde28uPW3RUGgDzD2ZDspMaEI7GQqrW_kqNHq3uTTG_hdx_o8e63dFZ0HXjWT-D09ziZd7_6A2p4NCf2aDNW_Ys8y0PkrVX1VyLHxqpT2j62-t9y5l-fd3AlRkgAe1WnmJVOUTKATHzYKK4MQKSW5bQ6fMV1Khb_HjVWzk-zpoZGBKNYkPEAGmQQv7Tmr9CZ_gD2aoCDibXCFBXi2k-TftSsnlzKbZkIF2IdSPHWpthXwj6Xf44Qnx1HRCp10MlrZmPI49TlNa2ZzHQUm2TEvFuqwzUMB7qIPWUEGdgSfIUpc00t3Yhx_2BaGbIzCp6Cp4O6BgTymoRSErj_Xg29Y4Xzkyen1nUc2zc0pWSd_h-iK5Jt9FmOupM4LJVSM_b04nIPARZFR7GxhRbzE5YX-EZhy2_yJIsH_Ajl0fPdZ2sDm3zne3OhdcBZ-Ea0BEqpemjOh1a717iX1b2LbLf5wkx2lmlVUhkO09FtwYqw-U26jrElqmfOl6vEDvG-GUJoYTHrSJ_nVFbaqt2eGqU6HGSyS7FQa0LVR5jANJPU2vi78s8y0FF75_CIMYHpjVYyQXkgEGNVmpm5wKBdagk4E9YSeSRey0gETkz1iVigAf80EL-c5nvMUfCvOktgDrjRo./download
    Resolving public.boxcloud.com (public.boxcloud.com)... 107.152.26.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|107.152.26.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 1775938697 (1.7G) [application/octet-stream]
    Saving to: ‘nuclei_900_S1_L001_R2_001.fastq.gz’
    
    nuclei_900_S1_L001_ 100%[===================>]   1.65G  25.4MB/s    in 66s     
    
    2020-01-16 15:53:50 (25.6 MB/s) - ‘nuclei_900_S1_L001_R2_001.fastq.gz’ saved [1775938697/1775938697]
    
    --2020-01-16 15:53:50--  https://caltech.box.com/shared/static/tlqdm0w3tvy8ogyktsz7ahggwurc6kkj.gz
    Resolving caltech.box.com (caltech.box.com)... 107.152.25.197, 107.152.24.197
    Connecting to caltech.box.com (caltech.box.com)|107.152.25.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/tlqdm0w3tvy8ogyktsz7ahggwurc6kkj.gz [following]
    --2020-01-16 15:53:50--  https://caltech.box.com/public/static/tlqdm0w3tvy8ogyktsz7ahggwurc6kkj.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/tlqdm0w3tvy8ogyktsz7ahggwurc6kkj.gz [following]
    --2020-01-16 15:53:50--  https://caltech.app.box.com/public/static/tlqdm0w3tvy8ogyktsz7ahggwurc6kkj.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 107.152.26.199, 107.152.27.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|107.152.26.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!UpSBIsWK59JE_USL_X0zT7AA4IRB7hbcmLNrL0FEo0Dt4ThSyOx7Q8_J2S8riTT6_4D3PABwaLWLqxyQBOh237ZHCJp7bjubDL7-JUopsW6dWCNgkwBeI5HQT_qya_pK7qd6kVvhP-YC_ssA7aueoRqMkPVnLRR2vEs1373v4ws6GQun7Cv4KhQsy3Y1G9fSL5a7jnoyhVAKY1T_jb-3bMlaEe0Lw5lOppJ1XvKLTasEVlmxJYBpsBUkGFMf2ReI5hBU4YCx-UlZiyKNFBs0-gwb6XjlbrY4VOSlvGwxYwC0AVRmabz3RPM0-grARJAlf2Nbmc_xSJjWKrlW8e79Ac4tSNG7IHEWX_O5wxLYAOfaLWMsZnnfGXAzE5a5IiYfhA7dq0egt8SBf2e5HraIwhkwO3_ud7hh0To2jcRJ-0c4VDk2eus9uYm92NDOSm0gxFSrSVIRpeK59y-V0OcjDDxjbVvHJxhKpEZET0G4q84NudifP8diW5sEKcvU4xTyUxF-oYJwv6HEvlXdtS5AQM-hLK-sjiJtdPU-8taL5f-xOkTIxnngzRRhlyDj8dy4eRDEgJlStPEdpUidCinmNsTCKvzqCuUcLoruGNUZawSSUwKqISaYkjNqqT22Tazth-Q2dnDCx4zo9btVH3Th9FZeff3JpsdtVRDkHLC-k8xcqGbAqG8bS7gDo8p3SeBDcCVzWibYqSwBYkGQ7WXs4Rgixwzd8SA2W-2A5pjZUnA03EEOMr6fCezBKA6u4VrpeKB_zfnDNwMPJCRouU7bkWPx6JmEczPfF7MR3EQW_HAEqVeM7L5g7YOxBq_TzgQz4kiwKENa97dNXPsFW96G8Wz2EkPzmfAHSvTM154ZNWB2H593eLU9EL1VskNZuUis-RFlsSNhXQkjcbbWJ-e48aKCU8tz67SqEaOs9jI-5BgblQp5sYl7KfquV3_VYUMqhB0fFrRW_5S88vHuje21yAl7TskomNDi0IPserkvsRECt7EbqIReDHeskpZY_YayU7I_CXjMXCIZ2UeqXGsI9Uw1kt0A55yGHTlu3nGb0D4lnXtVq9rA8WG483qATUmS23LYXbmbWa7_cJUWQclBC4fh-2C4h3oIwfOMrqgFUgwIW7FCScT8i1REUvtG_WaGv-79AziBs3R7loL-q1TR_HPMSkEbPT6lZWcXuF-orWadDZg9mcmxvIVNwwI2Kcbx5AnZa6BkuF1NRLyMc6kSC3A4E1ml8UXLmW7Ti_Hoc_WO6pTmnkkZd9ODWXM9D9eiBfbwfLExxgbws8Mt85EVqoV2fbSGTjlOeEsToWuhmhYfDSYeDZw-5-tjt1UnNoaD5D5z54CbnXZxX3uvx47Ninx_KnKvFpXN6EbPYo2PXYIO8il09FajIj3tNBXEIQZERszgeSpNGZkmcufJr8jprNP27w../download [following]
    --2020-01-16 15:53:51--  https://public.boxcloud.com/d/1/b1!UpSBIsWK59JE_USL_X0zT7AA4IRB7hbcmLNrL0FEo0Dt4ThSyOx7Q8_J2S8riTT6_4D3PABwaLWLqxyQBOh237ZHCJp7bjubDL7-JUopsW6dWCNgkwBeI5HQT_qya_pK7qd6kVvhP-YC_ssA7aueoRqMkPVnLRR2vEs1373v4ws6GQun7Cv4KhQsy3Y1G9fSL5a7jnoyhVAKY1T_jb-3bMlaEe0Lw5lOppJ1XvKLTasEVlmxJYBpsBUkGFMf2ReI5hBU4YCx-UlZiyKNFBs0-gwb6XjlbrY4VOSlvGwxYwC0AVRmabz3RPM0-grARJAlf2Nbmc_xSJjWKrlW8e79Ac4tSNG7IHEWX_O5wxLYAOfaLWMsZnnfGXAzE5a5IiYfhA7dq0egt8SBf2e5HraIwhkwO3_ud7hh0To2jcRJ-0c4VDk2eus9uYm92NDOSm0gxFSrSVIRpeK59y-V0OcjDDxjbVvHJxhKpEZET0G4q84NudifP8diW5sEKcvU4xTyUxF-oYJwv6HEvlXdtS5AQM-hLK-sjiJtdPU-8taL5f-xOkTIxnngzRRhlyDj8dy4eRDEgJlStPEdpUidCinmNsTCKvzqCuUcLoruGNUZawSSUwKqISaYkjNqqT22Tazth-Q2dnDCx4zo9btVH3Th9FZeff3JpsdtVRDkHLC-k8xcqGbAqG8bS7gDo8p3SeBDcCVzWibYqSwBYkGQ7WXs4Rgixwzd8SA2W-2A5pjZUnA03EEOMr6fCezBKA6u4VrpeKB_zfnDNwMPJCRouU7bkWPx6JmEczPfF7MR3EQW_HAEqVeM7L5g7YOxBq_TzgQz4kiwKENa97dNXPsFW96G8Wz2EkPzmfAHSvTM154ZNWB2H593eLU9EL1VskNZuUis-RFlsSNhXQkjcbbWJ-e48aKCU8tz67SqEaOs9jI-5BgblQp5sYl7KfquV3_VYUMqhB0fFrRW_5S88vHuje21yAl7TskomNDi0IPserkvsRECt7EbqIReDHeskpZY_YayU7I_CXjMXCIZ2UeqXGsI9Uw1kt0A55yGHTlu3nGb0D4lnXtVq9rA8WG483qATUmS23LYXbmbWa7_cJUWQclBC4fh-2C4h3oIwfOMrqgFUgwIW7FCScT8i1REUvtG_WaGv-79AziBs3R7loL-q1TR_HPMSkEbPT6lZWcXuF-orWadDZg9mcmxvIVNwwI2Kcbx5AnZa6BkuF1NRLyMc6kSC3A4E1ml8UXLmW7Ti_Hoc_WO6pTmnkkZd9ODWXM9D9eiBfbwfLExxgbws8Mt85EVqoV2fbSGTjlOeEsToWuhmhYfDSYeDZw-5-tjt1UnNoaD5D5z54CbnXZxX3uvx47Ninx_KnKvFpXN6EbPYo2PXYIO8il09FajIj3tNBXEIQZERszgeSpNGZkmcufJr8jprNP27w../download
    Resolving public.boxcloud.com (public.boxcloud.com)... 107.152.24.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|107.152.24.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 630554964 (601M) [application/octet-stream]
    Saving to: ‘nuclei_900_S1_L002_R1_001.fastq.gz’
    
    nuclei_900_S1_L002_ 100%[===================>] 601.34M  12.4MB/s    in 47s     
    
    2020-01-16 15:54:38 (12.8 MB/s) - ‘nuclei_900_S1_L002_R1_001.fastq.gz’ saved [630554964/630554964]
    
    --2020-01-16 15:54:39--  https://caltech.box.com/shared/static/gqrvkqllr9d7zq4e3yfrng9kgfbejowe.gz
    Resolving caltech.box.com (caltech.box.com)... 107.152.25.197, 107.152.24.197
    Connecting to caltech.box.com (caltech.box.com)|107.152.25.197|:443... connected.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: /public/static/gqrvkqllr9d7zq4e3yfrng9kgfbejowe.gz [following]
    --2020-01-16 15:54:39--  https://caltech.box.com/public/static/gqrvkqllr9d7zq4e3yfrng9kgfbejowe.gz
    Reusing existing connection to caltech.box.com:443.
    HTTP request sent, awaiting response... 301 Moved Permanently
    Location: https://caltech.app.box.com/public/static/gqrvkqllr9d7zq4e3yfrng9kgfbejowe.gz [following]
    --2020-01-16 15:54:39--  https://caltech.app.box.com/public/static/gqrvkqllr9d7zq4e3yfrng9kgfbejowe.gz
    Resolving caltech.app.box.com (caltech.app.box.com)... 107.152.26.199, 107.152.27.199
    Connecting to caltech.app.box.com (caltech.app.box.com)|107.152.26.199|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://public.boxcloud.com/d/1/b1!eSBRhEjnaGBNdIysLMxDJwbEGKBI6iXrGsTMIedrDwgQbgwehejXPFG-E1iLSBUDfYdCUl87Kl6-1c0h5LLPc6_S8Xv1GciJSIP4FQsSFnJI2yk304f87Z69O2sdOlHtOe4KFe5abT0R1reDd8eaXqT0QO8AkjttXlH1Nn4CJ7nlCXZWptCNrOyfmi5cRbLt-dnWHUyv3o1YzCnoYAvvW-lyXcQloP25G80fpfzHDry8e5e4WlDZag9AOrtrmc8t3QLH8875LU3-AM-VT3Cyyca5LFcZJ42ESv62w6fnQHongPS9Rw-1ULp3nytnmQLgWF4E4kxYst7qUYV5Jmfg1UzboCi2ld0nBT4FWhpm8r37k51ZDNYc338UWJizAsjjHtO4e0vrAukWOgOPwyof7_jVx-g6wj0f7BXu3LNFMbbl9c6DkAsTGP9k-dfb9uELb1iLiYQQTSsxw9FXKJoxepRZJVPKnfUV4QOy449S8KJlCqe_gopt0u6V7E7yya9GhPIZtgtkgVLGSGCgaGsSETTy4Ab3hnOIWVhwocPEr_AST5IcC1WI71q22n6gY80eMXa61Pu0bzCI-gIoKGIGrnX8D1SaYfUK5EXcZ2dfhVrU85Dj_FuOK4_QFsFVRzbJfzowewSFOCIKt6kQDxvTcz2gIC5puHg8DJL0thNSY5YonBMHbbuikl7nlr2S2NT24X9ii7i-cMX87abVNpu1gir8D0qtRzrdUNc5-1UMMmXfyrGT5frLxPT0wZZyfvVHSvoIf22EDBMD4_h4hlt_MGXExuPhsroeIaBhMBKnEJkp5mOsFJsKfGh8ZIbAwDu0Q3GPKFdIqY1LWuwL2qOQDqjG2KM2FDk_bh-qNy9xVQ8swIYXMlEB-9mC5wXtsOAHJzb9f2q9P3dqn0Jc-vYINs-onNOYwGqHHfGyupQJ3GqAdvVWcuhy-ARGj4XLVied2vhgC8E1YxlQ9jlZMszuwt-3_fLKZ_xwlz4VTPLtP7hNYc9OIEKg03yViIiyBFNllviCVHs1ej0Xicmn5S5Z_R8toCOHUkS4hiz560zZbCwZk5BmIqYdne3Im5cOeeDN_FNVxxL9njXlj7FZm7p1q_sKOxw0oUHGp54URZEuFWP6o4HgB-Bo0ZzWyTjgiJWHzww0UJ8uUiphlInz4FLc72k7lx2En04_SHlP2VrXF7nAdkaSQOynlerPO03bPfzmDNE1ZDTq_KIBD9W_9WD36nmOHJKY34kc-SWR7s-g8494ASxsIvubNonM1wEivcUJt3pazsE3IMZNuo2_cfh0DcARwuoCAunGR18wBgSh85XdpUwI5E6WXylPdyphNfNiso2aX3MIk7018KlJJ-vzrDG7U3bArSmgmtUnd1vzVKw0a1a4QHCyeKgJR0FyaP0vAjj5KQNpZAcfd0lVqfdPQOaPa5AY/download [following]
    --2020-01-16 15:54:40--  https://public.boxcloud.com/d/1/b1!eSBRhEjnaGBNdIysLMxDJwbEGKBI6iXrGsTMIedrDwgQbgwehejXPFG-E1iLSBUDfYdCUl87Kl6-1c0h5LLPc6_S8Xv1GciJSIP4FQsSFnJI2yk304f87Z69O2sdOlHtOe4KFe5abT0R1reDd8eaXqT0QO8AkjttXlH1Nn4CJ7nlCXZWptCNrOyfmi5cRbLt-dnWHUyv3o1YzCnoYAvvW-lyXcQloP25G80fpfzHDry8e5e4WlDZag9AOrtrmc8t3QLH8875LU3-AM-VT3Cyyca5LFcZJ42ESv62w6fnQHongPS9Rw-1ULp3nytnmQLgWF4E4kxYst7qUYV5Jmfg1UzboCi2ld0nBT4FWhpm8r37k51ZDNYc338UWJizAsjjHtO4e0vrAukWOgOPwyof7_jVx-g6wj0f7BXu3LNFMbbl9c6DkAsTGP9k-dfb9uELb1iLiYQQTSsxw9FXKJoxepRZJVPKnfUV4QOy449S8KJlCqe_gopt0u6V7E7yya9GhPIZtgtkgVLGSGCgaGsSETTy4Ab3hnOIWVhwocPEr_AST5IcC1WI71q22n6gY80eMXa61Pu0bzCI-gIoKGIGrnX8D1SaYfUK5EXcZ2dfhVrU85Dj_FuOK4_QFsFVRzbJfzowewSFOCIKt6kQDxvTcz2gIC5puHg8DJL0thNSY5YonBMHbbuikl7nlr2S2NT24X9ii7i-cMX87abVNpu1gir8D0qtRzrdUNc5-1UMMmXfyrGT5frLxPT0wZZyfvVHSvoIf22EDBMD4_h4hlt_MGXExuPhsroeIaBhMBKnEJkp5mOsFJsKfGh8ZIbAwDu0Q3GPKFdIqY1LWuwL2qOQDqjG2KM2FDk_bh-qNy9xVQ8swIYXMlEB-9mC5wXtsOAHJzb9f2q9P3dqn0Jc-vYINs-onNOYwGqHHfGyupQJ3GqAdvVWcuhy-ARGj4XLVied2vhgC8E1YxlQ9jlZMszuwt-3_fLKZ_xwlz4VTPLtP7hNYc9OIEKg03yViIiyBFNllviCVHs1ej0Xicmn5S5Z_R8toCOHUkS4hiz560zZbCwZk5BmIqYdne3Im5cOeeDN_FNVxxL9njXlj7FZm7p1q_sKOxw0oUHGp54URZEuFWP6o4HgB-Bo0ZzWyTjgiJWHzww0UJ8uUiphlInz4FLc72k7lx2En04_SHlP2VrXF7nAdkaSQOynlerPO03bPfzmDNE1ZDTq_KIBD9W_9WD36nmOHJKY34kc-SWR7s-g8494ASxsIvubNonM1wEivcUJt3pazsE3IMZNuo2_cfh0DcARwuoCAunGR18wBgSh85XdpUwI5E6WXylPdyphNfNiso2aX3MIk7018KlJJ-vzrDG7U3bArSmgmtUnd1vzVKw0a1a4QHCyeKgJR0FyaP0vAjj5KQNpZAcfd0lVqfdPQOaPa5AY/download
    Resolving public.boxcloud.com (public.boxcloud.com)... 107.152.24.200
    Connecting to public.boxcloud.com (public.boxcloud.com)|107.152.24.200|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 1784669831 (1.7G) [application/octet-stream]
    Saving to: ‘nuclei_900_S1_L002_R2_001.fastq.gz’
    
    nuclei_900_S1_L002_ 100%[===================>]   1.66G  17.5MB/s    in 95s     
    
    2020-01-16 15:56:15 (17.9 MB/s) - ‘nuclei_900_S1_L002_R2_001.fastq.gz’ saved [1784669831/1784669831]
    
    CPU times: user 2.09 s, sys: 407 ms, total: 2.49 s
    Wall time: 4min 9s


Then, we verify the integrity of the files we downloaded to make sure they were not corrupted during the download.


```
!md5sum -c checksums.txt --ignore-missing
```

    nuclei_900_S1_L001_R1_001.fastq.gz: OK
    nuclei_900_S1_L001_R2_001.fastq.gz: OK
    nuclei_900_S1_L002_R1_001.fastq.gz: OK
    nuclei_900_S1_L002_R2_001.fastq.gz: OK


### Install `kb`

Install `kb` for running the **kallisto | bustools** workflow.


```
!pip install git+https://github.com/pachterlab/kb_python@count-kite
```

    Collecting git+https://github.com/pachterlab/kb_python@count-kite
      Cloning https://github.com/pachterlab/kb_python (to revision count-kite) to /tmp/pip-req-build-00ne_j5p
      Running command git clone -q https://github.com/pachterlab/kb_python /tmp/pip-req-build-00ne_j5p
      Running command git checkout -b count-kite --track origin/count-kite
      Switched to a new branch 'count-kite'
      Branch 'count-kite' set up to track remote branch 'count-kite' from 'origin'.
    Collecting anndata>=0.6.22.post1
    [?25l  Downloading https://files.pythonhosted.org/packages/2b/72/87196c15f68d9865c31a43a10cf7c50bcbcedd5607d09f9aada0b3963103/anndata-0.6.22.post1-py3-none-any.whl (47kB)
    [K     |████████████████████████████████| 51kB 1.8MB/s 
    [?25hCollecting loompy>=3.0.6
    [?25l  Downloading https://files.pythonhosted.org/packages/36/52/74ed37ae5988522fbf87b856c67c4f80700e6452410b4cd80498c5f416f9/loompy-3.0.6.tar.gz (41kB)
    [K     |████████████████████████████████| 51kB 5.1MB/s 
    [?25hRequirement already satisfied: requests>=2.19.0 in /usr/local/lib/python3.6/dist-packages (from kb-python==0.24.4) (2.21.0)
    Collecting tqdm>=4.39.0
    [?25l  Downloading https://files.pythonhosted.org/packages/72/c9/7fc20feac72e79032a7c8138fd0d395dc6d8812b5b9edf53c3afd0b31017/tqdm-4.41.1-py2.py3-none-any.whl (56kB)
    [K     |████████████████████████████████| 61kB 3.9MB/s 
    [?25hRequirement already satisfied: pandas>=0.23.0 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (0.25.3)
    Requirement already satisfied: numpy~=1.14 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (1.17.5)
    Requirement already satisfied: h5py in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (2.8.0)
    Requirement already satisfied: natsort in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (5.5.0)
    Requirement already satisfied: scipy~=1.0 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (1.4.1)
    Requirement already satisfied: setuptools in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python==0.24.4) (42.0.2)
    Requirement already satisfied: numba in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python==0.24.4) (0.47.0)
    Requirement already satisfied: click in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python==0.24.4) (7.0)
    Collecting numpy-groupies
    [?25l  Downloading https://files.pythonhosted.org/packages/57/ae/18217b57ba3e4bb8a44ecbfc161ed065f6d1b90c75d404bd6ba8d6f024e2/numpy_groupies-0.9.10.tar.gz (43kB)
    [K     |████████████████████████████████| 51kB 4.7MB/s 
    [?25hRequirement already satisfied: urllib3<1.25,>=1.21.1 in /usr/local/lib/python3.6/dist-packages (from requests>=2.19.0->kb-python==0.24.4) (1.24.3)
    Requirement already satisfied: certifi>=2017.4.17 in /usr/local/lib/python3.6/dist-packages (from requests>=2.19.0->kb-python==0.24.4) (2019.11.28)
    Requirement already satisfied: idna<2.9,>=2.5 in /usr/local/lib/python3.6/dist-packages (from requests>=2.19.0->kb-python==0.24.4) (2.8)
    Requirement already satisfied: chardet<3.1.0,>=3.0.2 in /usr/local/lib/python3.6/dist-packages (from requests>=2.19.0->kb-python==0.24.4) (3.0.4)
    Requirement already satisfied: python-dateutil>=2.6.1 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python==0.24.4) (2.6.1)
    Requirement already satisfied: pytz>=2017.2 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python==0.24.4) (2018.9)
    Requirement already satisfied: six in /usr/local/lib/python3.6/dist-packages (from h5py->anndata>=0.6.22.post1->kb-python==0.24.4) (1.12.0)
    Requirement already satisfied: llvmlite>=0.31.0dev0 in /usr/local/lib/python3.6/dist-packages (from numba->loompy>=3.0.6->kb-python==0.24.4) (0.31.0)
    Building wheels for collected packages: kb-python, loompy, numpy-groupies
      Building wheel for kb-python (setup.py) ... [?25l[?25hdone
      Created wheel for kb-python: filename=kb_python-0.24.4-cp36-none-any.whl size=80991434 sha256=46c80d686746b8cc5f5ea19481297c38fcb7c928c8f472d75d4261ce74266483
      Stored in directory: /tmp/pip-ephem-wheel-cache-6ek0s5dl/wheels/8e/56/56/c89223de74af26792675e82f4bb5223e7cf0d653a33038e34c
      Building wheel for loompy (setup.py) ... [?25l[?25hdone
      Created wheel for loompy: filename=loompy-3.0.6-cp36-none-any.whl size=47896 sha256=5aa0e8271e5ecd4bdfc5ad2636af749cb6b32c6dc5789e4cca70337993c1f066
      Stored in directory: /root/.cache/pip/wheels/f9/a4/90/5a98ad83419732b0fba533b81a2a52ba3dbe230a936ca4cdc9
      Building wheel for numpy-groupies (setup.py) ... [?25l[?25hdone
      Created wheel for numpy-groupies: filename=numpy_groupies-0+unknown-cp36-none-any.whl size=28044 sha256=7f56b834c065c60c4db18f99a2d651074b23e468028919e5ff9422f774ebc9cf
      Stored in directory: /root/.cache/pip/wheels/30/ac/83/64d5f9293aeaec63f9539142fc629a41af064cae1b3d8d94aa
    Successfully built kb-python loompy numpy-groupies
    Installing collected packages: anndata, numpy-groupies, loompy, tqdm, kb-python
      Found existing installation: tqdm 4.28.1
        Uninstalling tqdm-4.28.1:
          Successfully uninstalled tqdm-4.28.1
    Successfully installed anndata-0.6.22.post1 kb-python-0.24.4 loompy-3.0.6 numpy-groupies-0+unknown tqdm-4.41.1




### Download mouse reference files and build the index

We build a mouse cDNA and intron index from the mouse genome and annotations provided by Ensembl.


```
%%time
!wget ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
!wget ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz
```

    --2020-01-16 15:56:57--  ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
               => ‘Mus_musculus.GRCm38.dna.primary_assembly.fa.gz’
    Resolving ftp.ensembl.org (ftp.ensembl.org)... 193.62.193.8
    Connecting to ftp.ensembl.org (ftp.ensembl.org)|193.62.193.8|:21... connected.
    Logging in as anonymous ... Logged in!
    ==> SYST ... done.    ==> PWD ... done.
    ==> TYPE I ... done.  ==> CWD (1) /pub/release-98/fasta/mus_musculus/dna ... done.
    ==> SIZE Mus_musculus.GRCm38.dna.primary_assembly.fa.gz ... 805984352
    ==> PASV ... done.    ==> RETR Mus_musculus.GRCm38.dna.primary_assembly.fa.gz ... done.
    Length: 805984352 (769M) (unauthoritative)
    
    Mus_musculus.GRCm38 100%[===================>] 768.65M  23.0MB/s    in 36s     
    
    2020-01-16 15:57:35 (21.1 MB/s) - ‘Mus_musculus.GRCm38.dna.primary_assembly.fa.gz’ saved [805984352]
    
    --2020-01-16 15:57:36--  ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz
               => ‘Mus_musculus.GRCm38.98.gtf.gz’
    Resolving ftp.ensembl.org (ftp.ensembl.org)... 193.62.193.8
    Connecting to ftp.ensembl.org (ftp.ensembl.org)|193.62.193.8|:21... connected.
    Logging in as anonymous ... Logged in!
    ==> SYST ... done.    ==> PWD ... done.
    ==> TYPE I ... done.  ==> CWD (1) /pub/release-98/gtf/mus_musculus ... done.
    ==> SIZE Mus_musculus.GRCm38.98.gtf.gz ... 30256597
    ==> PASV ... done.    ==> RETR Mus_musculus.GRCm38.98.gtf.gz ... done.
    Length: 30256597 (29M) (unauthoritative)
    
    Mus_musculus.GRCm38 100%[===================>]  28.85M  14.0MB/s    in 2.1s    
    
    2020-01-16 15:57:40 (14.0 MB/s) - ‘Mus_musculus.GRCm38.98.gtf.gz’ saved [30256597]
    
    CPU times: user 371 ms, sys: 76.5 ms, total: 447 ms
    Wall time: 43.2 s


Notice we use the `-n` option to split the index into multiple files. This helps reduce the maximum memory usage so that Google Colab doesn't run out of memory.


```
%%time
!kb ref -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa \
-c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno -n 8 \
Mus_musculus.GRCm38.dna.primary_assembly.fa.gz Mus_musculus.GRCm38.98.gtf.gz
```

    [2020-01-16 15:57:41,860]    INFO Preparing Mus_musculus.GRCm38.dna.primary_assembly.fa.gz, Mus_musculus.GRCm38.98.gtf.gz
    [2020-01-16 15:57:41,860]    INFO Decompressing Mus_musculus.GRCm38.dna.primary_assembly.fa.gz to tmp
    [2020-01-16 15:58:08,556]    INFO Sorting tmp/Mus_musculus.GRCm38.dna.primary_assembly.fa to /content/tmp/tmptuipvatl
    [2020-01-16 16:05:44,930]    INFO Decompressing Mus_musculus.GRCm38.98.gtf.gz to tmp
    [2020-01-16 16:05:48,839]    INFO Sorting tmp/Mus_musculus.GRCm38.98.gtf to /content/tmp/tmpqc7d85yg
    [2020-01-16 16:06:50,352]    INFO Splitting genome tmp/Mus_musculus.GRCm38.dna.primary_assembly.fa into cDNA at /content/tmp/tmpzpi_mdzg
    [2020-01-16 16:06:50,352] WARNING The following chromosomes were found in the FASTA but doens't have any "transcript" features in the GTF: GL456383.1, GL456367.1, GL456390.1, GL456389.1, GL456379.1, GL456387.1, GL456370.1, GL456394.1, GL456359.1, GL456392.1, GL456378.1, GL456213.1, GL456396.1, JH584300.1, GL456360.1, GL456393.1, JH584301.1, GL456366.1, JH584302.1, GL456382.1, GL456368.1. No sequences will be generated for these chromosomes.
    [2020-01-16 16:07:59,570]    INFO Wrote 142446 cDNA transcripts
    [2020-01-16 16:07:59,575]    INFO Creating cDNA transcripts-to-capture at /content/tmp/tmpgun332i9
    [2020-01-16 16:08:00,603]    INFO Splitting genome into introns at /content/tmp/tmp4cl5pq5u
    [2020-01-16 16:12:49,090]    INFO Wrote 647972 intron sequences
    [2020-01-16 16:12:49,095]    INFO Creating intron transcripts-to-capture at /content/tmp/tmpb7kwg4ry
    [2020-01-16 16:13:55,504]    INFO Concatenating 1 cDNA FASTAs to cdna.fa
    [2020-01-16 16:14:00,221]    INFO Concatenating 1 cDNA transcripts-to-captures to cdna_t2c.txt
    [2020-01-16 16:14:00,320]    INFO Concatenating 1 intron FASTAs to intron.fa
    [2020-01-16 16:14:50,492]    INFO Concatenating 1 intron transcripts-to-captures to intron_t2c.txt
    [2020-01-16 16:14:50,803]    INFO Concatenating cDNA and intron FASTAs to /content/tmp/tmp6ssxsc_j
    [2020-01-16 16:16:10,955]    INFO Creating transcript-to-gene mapping at t2g.txt
    [2020-01-16 16:17:18,925]    INFO Splitting /content/tmp/tmp6ssxsc_j into 8 parts
    [2020-01-16 16:18:13,532]    INFO Indexing /content/tmp/tmpzw65y9wa to index.idx.0
    [2020-01-16 16:35:27,873]    INFO Indexing /content/tmp/tmp4h1a02ef to index.idx.1
    [2020-01-16 16:49:46,788]    INFO Indexing /content/tmp/tmptw5rlyf0 to index.idx.2
    [2020-01-16 17:04:53,914]    INFO Indexing /content/tmp/tmpavi28g7n to index.idx.3
    [2020-01-16 17:19:13,835]    INFO Indexing /content/tmp/tmpofcnmlri to index.idx.4
    [2020-01-16 17:33:31,933]    INFO Indexing /content/tmp/tmpq_weqj_3 to index.idx.5
    [2020-01-16 17:48:06,886]    INFO Indexing /content/tmp/tmpp7xih6mo to index.idx.6
    [2020-01-16 18:03:09,010]    INFO Indexing /content/tmp/tmpyj0138gw to index.idx.7
    CPU times: user 36.7 s, sys: 6.37 s, total: 43.1 s
    Wall time: 2h 20min 44s


### Generate RNA count matrices

The following command will generate an RNA count matrix of cells (rows) by genes (columns) in H5AD format, which is a binary format used to store [Anndata](https://anndata.readthedocs.io/en/stable/) objects. Notice we are providing the index and transcript-to-gene mapping we downloaded in the previous step to the `-i` and `-g` arguments respectively, as well as the transcripts-to-capture lists to the `-c1` and `-c2` arguments. Also, these reads were generated with the 10x Genomics Chromium Single Cell v2 Chemistry, hence the `-x 10xv2` argument. To view other supported technologies, run `kb --list`. Note the `--workflow nucleus` to indicate this is a single nucleus experiment.

__Note:__ If you would like a Loom file instead, replace the `--h5ad` flag with `--loom`. If you want to use the raw matrix output by `kb` instead of their H5AD or Loom converted files, omit these flags.


```
%%time
!kb count -i index.idx.0,index.idx.1,index.idx.2,index.idx.3,index.idx.4,index.idx.5,index.idx.6,index.idx.7 \
-g t2g.txt -c1 cdna_t2c.txt -c2 intron_t2c.txt -x 10xv2 -o output -t 2 --workflow nucleus --h5ad \
nuclei_900_S1_L001_R1_001.fastq.gz nuclei_900_S1_L001_R2_001.fastq.gz \
nuclei_900_S1_L002_R1_001.fastq.gz nuclei_900_S1_L002_R2_001.fastq.gz
```

    [2020-01-16 21:26:55,274]    INFO Generating BUS file using 8 indices
    [2020-01-16 21:26:55,274]    INFO Generating BUS file to output/tmp/bus_part0 from
    [2020-01-16 21:26:55,274]    INFO         nuclei_900_S1_L001_R1_001.fastq.gz
    [2020-01-16 21:26:55,274]    INFO         nuclei_900_S1_L001_R2_001.fastq.gz
    [2020-01-16 21:26:55,274]    INFO         nuclei_900_S1_L002_R1_001.fastq.gz
    [2020-01-16 21:26:55,274]    INFO         nuclei_900_S1_L002_R2_001.fastq.gz
    [2020-01-16 21:26:55,275]    INFO Using index index.idx.0
    [2020-01-16 21:38:46,592]    INFO Generating BUS file to output/tmp/bus_part1 from
    [2020-01-16 21:38:46,592]    INFO         nuclei_900_S1_L001_R1_001.fastq.gz
    [2020-01-16 21:38:46,592]    INFO         nuclei_900_S1_L001_R2_001.fastq.gz
    [2020-01-16 21:38:46,593]    INFO         nuclei_900_S1_L002_R1_001.fastq.gz
    [2020-01-16 21:38:46,593]    INFO         nuclei_900_S1_L002_R2_001.fastq.gz
    [2020-01-16 21:38:46,593]    INFO Using index index.idx.1
    [2020-01-16 21:54:34,715]    INFO Generating BUS file to output/tmp/bus_part2 from
    [2020-01-16 21:54:34,716]    INFO         nuclei_900_S1_L001_R1_001.fastq.gz
    [2020-01-16 21:54:34,716]    INFO         nuclei_900_S1_L001_R2_001.fastq.gz
    [2020-01-16 21:54:34,716]    INFO         nuclei_900_S1_L002_R1_001.fastq.gz
    [2020-01-16 21:54:34,716]    INFO         nuclei_900_S1_L002_R2_001.fastq.gz
    [2020-01-16 21:54:34,716]    INFO Using index index.idx.2
    [2020-01-16 22:10:18,962]    INFO Generating BUS file to output/tmp/bus_part3 from
    [2020-01-16 22:10:18,963]    INFO         nuclei_900_S1_L001_R1_001.fastq.gz
    [2020-01-16 22:10:18,963]    INFO         nuclei_900_S1_L001_R2_001.fastq.gz
    [2020-01-16 22:10:18,963]    INFO         nuclei_900_S1_L002_R1_001.fastq.gz
    [2020-01-16 22:10:18,963]    INFO         nuclei_900_S1_L002_R2_001.fastq.gz
    [2020-01-16 22:10:18,963]    INFO Using index index.idx.3
    [2020-01-16 22:26:04,028]    INFO Generating BUS file to output/tmp/bus_part4 from
    [2020-01-16 22:26:04,028]    INFO         nuclei_900_S1_L001_R1_001.fastq.gz
    [2020-01-16 22:26:04,029]    INFO         nuclei_900_S1_L001_R2_001.fastq.gz
    [2020-01-16 22:26:04,029]    INFO         nuclei_900_S1_L002_R1_001.fastq.gz
    [2020-01-16 22:26:04,029]    INFO         nuclei_900_S1_L002_R2_001.fastq.gz
    [2020-01-16 22:26:04,029]    INFO Using index index.idx.4
    [2020-01-16 22:41:47,809]    INFO Generating BUS file to output/tmp/bus_part5 from
    [2020-01-16 22:41:47,809]    INFO         nuclei_900_S1_L001_R1_001.fastq.gz
    [2020-01-16 22:41:47,809]    INFO         nuclei_900_S1_L001_R2_001.fastq.gz
    [2020-01-16 22:41:47,809]    INFO         nuclei_900_S1_L002_R1_001.fastq.gz
    [2020-01-16 22:41:47,809]    INFO         nuclei_900_S1_L002_R2_001.fastq.gz
    [2020-01-16 22:41:47,809]    INFO Using index index.idx.5
    [2020-01-16 22:57:47,130]    INFO Generating BUS file to output/tmp/bus_part6 from
    [2020-01-16 22:57:47,130]    INFO         nuclei_900_S1_L001_R1_001.fastq.gz
    [2020-01-16 22:57:47,130]    INFO         nuclei_900_S1_L001_R2_001.fastq.gz
    [2020-01-16 22:57:47,130]    INFO         nuclei_900_S1_L002_R1_001.fastq.gz
    [2020-01-16 22:57:47,130]    INFO         nuclei_900_S1_L002_R2_001.fastq.gz
    [2020-01-16 22:57:47,130]    INFO Using index index.idx.6
    [2020-01-16 23:13:46,707]    INFO Generating BUS file to output/tmp/bus_part7 from
    [2020-01-16 23:13:46,708]    INFO         nuclei_900_S1_L001_R1_001.fastq.gz
    [2020-01-16 23:13:46,708]    INFO         nuclei_900_S1_L001_R2_001.fastq.gz
    [2020-01-16 23:13:46,708]    INFO         nuclei_900_S1_L002_R1_001.fastq.gz
    [2020-01-16 23:13:46,708]    INFO         nuclei_900_S1_L002_R2_001.fastq.gz
    [2020-01-16 23:13:46,708]    INFO Using index index.idx.7
    [2020-01-16 23:29:45,121]    INFO Merging BUS records to output from
    [2020-01-16 23:29:45,121]    INFO         output/tmp/bus_part0
    [2020-01-16 23:29:45,121]    INFO         output/tmp/bus_part1
    [2020-01-16 23:29:45,121]    INFO         output/tmp/bus_part2
    [2020-01-16 23:29:45,121]    INFO         output/tmp/bus_part3
    [2020-01-16 23:29:45,121]    INFO         output/tmp/bus_part4
    [2020-01-16 23:29:45,121]    INFO         output/tmp/bus_part5
    [2020-01-16 23:29:45,121]    INFO         output/tmp/bus_part6
    [2020-01-16 23:29:45,121]    INFO         output/tmp/bus_part7
    [2020-01-16 23:33:12,815]    INFO Sorting BUS file output/output.bus to output/tmp/output.s.bus
    [2020-01-16 23:35:04,359]    INFO Whitelist not provided
    [2020-01-16 23:35:04,359]    INFO Copying pre-packaged 10XV2 whitelist to output
    [2020-01-16 23:35:04,517]    INFO Inspecting BUS file output/tmp/output.s.bus
    [2020-01-16 23:36:04,780]    INFO Correcting BUS records in output/tmp/output.s.bus to output/tmp/output.s.c.bus with whitelist output/10xv2_whitelist.txt
    [2020-01-16 23:37:08,241]    INFO Sorting BUS file output/tmp/output.s.c.bus to output/output.unfiltered.bus
    [2020-01-16 23:38:44,109]    INFO Capturing records from BUS file output/output.unfiltered.bus to output/tmp/spliced.bus with capture list intron_t2c.txt
    [2020-01-16 23:40:15,469]    INFO Sorting BUS file output/tmp/spliced.bus to output/spliced.unfiltered.bus
    [2020-01-16 23:40:31,065]    INFO Generating count matrix output/counts_unfiltered/spliced from BUS file output/spliced.unfiltered.bus
    [2020-01-16 23:41:40,242]    INFO Capturing records from BUS file output/output.unfiltered.bus to output/tmp/unspliced.bus with capture list cdna_t2c.txt
    [2020-01-16 23:43:55,755]    INFO Sorting BUS file output/tmp/unspliced.bus to output/unspliced.unfiltered.bus
    [2020-01-16 23:44:28,286]    INFO Generating count matrix output/counts_unfiltered/unspliced from BUS file output/unspliced.unfiltered.bus
    [2020-01-16 23:45:40,436]    INFO Reading matrix output/counts_unfiltered/spliced.mtx
    [2020-01-16 23:45:51,189]    INFO Reading matrix output/counts_unfiltered/unspliced.mtx
    [2020-01-16 23:45:57,929]    INFO Combining matrices
    [2020-01-16 23:45:58,596]    INFO Writing matrices to h5ad output/counts_unfiltered/adata.h5ad
    CPU times: user 38.4 s, sys: 5.08 s, total: 43.5 s
    Wall time: 2h 19min 8s


## Load the anndata into python

`kb` automatically sums the counts of genes for barcodes common to both spliced and unspliced matrices as an Anndata object in H5AD format (because `--h5ad` was used).


```
import anndata
```


```
adata = anndata.read('output/counts_unfiltered/adata.h5ad')
```


```
adata
```




    AnnData object with n_obs × n_vars = 175538 × 55421 




```
adata.obs
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
      <th>AAACCTGAGAAAGTGG</th>
    </tr>
    <tr>
      <th>AAACCTGAGAAGCCCA</th>
    </tr>
    <tr>
      <th>AAACCTGAGAATTCCC</th>
    </tr>
    <tr>
      <th>AAACCTGAGACACGAC</th>
    </tr>
    <tr>
      <th>AAACCTGAGACCCACC</th>
    </tr>
    <tr>
      <th>...</th>
    </tr>
    <tr>
      <th>TTTGTCATCTTACCGC</th>
    </tr>
    <tr>
      <th>TTTGTCATCTTAGCCC</th>
    </tr>
    <tr>
      <th>TTTGTCATCTTCAACT</th>
    </tr>
    <tr>
      <th>TTTGTCATCTTGAGGT</th>
    </tr>
    <tr>
      <th>TTTGTCATCTTTAGTC</th>
    </tr>
  </tbody>
</table>
<p>175538 rows × 0 columns</p>
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
      <th>ENSMUSG00000026535.9</th>
    </tr>
    <tr>
      <th>ENSMUSG00000026315.13</th>
    </tr>
    <tr>
      <th>ENSMUSG00000000817.10</th>
    </tr>
    <tr>
      <th>ENSMUSG00000063558.4</th>
    </tr>
    <tr>
      <th>ENSMUSG00000001138.13</th>
    </tr>
    <tr>
      <th>...</th>
    </tr>
    <tr>
      <th>ENSMUSG00000118447.1</th>
    </tr>
    <tr>
      <th>ENSMUSG00000118422.1</th>
    </tr>
    <tr>
      <th>ENSMUSG00000118472.1</th>
    </tr>
    <tr>
      <th>ENSMUSG00000118404.1</th>
    </tr>
    <tr>
      <th>ENSMUSG00000118417.1</th>
    </tr>
  </tbody>
</table>
<p>55421 rows × 0 columns</p>
</div>




```

```
