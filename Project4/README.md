#자유 프로젝트 
연구주제: Virus-infected human RNA-seq 데이터의 전사체 분석

Week1 : 
1. 논문을 읽고 Virus-infected human의 전사체 데이터를 받았습니다. 
  1)논문: “Deregulation of HDAC5 by Viral Interferon Regulatory Factor 3 Plays an Essential Role in Kaposi’s Sarcoma-Associated Herpesvirus-Induced Lymphangiogenesis”
2. 전사체 데이터의 전처리 파이프라인을 제작하여 count table형태로 만들었습니다.
   1)Trimmomatic-HISAT2-Samtools-featureCounts로 이어지는 파이프라인 제작

Week2: 
1. Raw data, TPM, TMM의 상관관계를 plotting
2. EdgeR를 이용하여 LEC(Lymphatic Epithelial Cell)과 vIRF3-overexpressed LEC 사이의 Differential Expressed Gene분석
3. Volcano plot을 plotting

Week3:
1. g:Profiler를 이용하여 logFC>4 & FDR<0.05인 DEG들에 대해 Gene Set Enrichment Analysis
2. Gene set enrichment analysis 결과 이해하기
