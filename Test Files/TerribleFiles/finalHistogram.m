function [  ] = finalHistogram(  )

bins = 40;

names=['Allen','An','Ayala','Bae','Barseghian','Belmonte','Bhakta','Brown','Burnham','Carrasco','Chen','Choe'...
'Conroy','Cosyn','Cuyno','Desopo','Doctor','Driscoll','Elia','Faulstick','Garcia','Gerlek','Ha'...
'Hadimulia','Hammock','Haughton','Hawkins','Horrocks','Hsu','Huang','Jiang','Kansal','Kolwalkar','Krief'...
'Kruger','Lau','Lazernik','Li','Litre','Liu','Liu','Lopez','Low','Manson','Matalon','Mesri','Nguyen'...
'Nguyen','Nguyen','Northway','Ortiz-martinez','Pacheco','Palaniswamy','Peng','Reyes','Reyes'...
'Sinnis','Smith','Sun','Thach','Tran','Tran','Tsai','Valenze','Van-linge','Velazquez-olivera','Villarreal'...
'Wang','Wang','Whitcomb','Williams','Xu','Young','Zhang','Zhu'];

arr=[66.30676
79.306
48.11054
68.539
58.59909
5.72094
70.19583
68.40743
47.57811
51.93943
85.98447
68.80643
61.71913
40.61116
52.04921
83.06573
64.05516
60.07122
79.93267
57.75347
65.41602
48.03231
65.72187
56.52412
93.44061
61.19171
47.97553
84.56221
52.60733
66.01553
69.70052
97.85638
49.1934
50.39276
58.32779
63.01208
35.44467
74.50691
70.09021
72.43229
81.23085
53.04799
60.39508
78.46462
70.72817
66.19058
51.395
61.22085
65.10936
68.6992
54.68344
69.55268
69.47908
72.30022
63.57849
25.95262
79.77526
49.73569
68.03553
48.84318
63.53748
81.69409
68.73916
56.12001
61.44043
79.3239
68.43975
82.81609
85.38676
51.00494
62.82723
79.05631
54.43603
81.17766
85.53782];

hist(arr,bins)

Average=mean(arr)
Standard_Dev=std(arr)

sum(arr>90)

end

