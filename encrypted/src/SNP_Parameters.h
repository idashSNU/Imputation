#ifndef SNPPARAMETERS_H_
#define SNPPARAMETERS_H_

class SNP_Parameters {
public:
    long list_x1k_start[101] = { 0, 137, 253, 283, 347, 415, 563, 646, 785, 894, 995, 1126, 1228, 1300, 1380, 1538, 1608, 1737, 1784, 1835, 1915, 2098, 2171, 2263, 2328, 2405, 2474, 2647, 2782, 2831, 2897, 3061, 3174, 3353, 3400, 3618, 3848, 4023, 4056, 4123, 4202, 4256, 4312, 4443, 4508, 4574, 4717, 4798, 4898, 5077, 5178, 5288, 5398, 5466, 5540, 5663, 5778, 5839, 5867, 5927, 5948, 6022, 6116, 6335, 6453, 6498, 6580, 6658, 6770, 6892, 6970, 7051, 7124, 7162, 7281, 7323, 7469, 7598, 7736, 7775, 7861, 7958, 8053, 8113, 8342, 8464, 8529, 8708, 8847, 8898, 8951, 9052, 9151, 9194, 9242, 9254, 9341, 9454, 9588, 9630, 9740 };

    long list_x1k_end[101] = { 136, 252, 282, 346, 414, 562, 645, 785, 893, 994, 1125, 1227, 1299, 1379, 1537, 1607, 1736, 1783, 1834, 1914, 2097, 2170, 2262, 2327, 2404, 2473, 2647, 2781, 2830, 2896, 3060, 3173, 3352, 3399, 3617, 3847, 4022, 4055, 4122, 4201, 4255, 4311, 4442, 4507, 4573, 4716, 4797, 4897, 5076, 5177, 5287, 5397, 5465, 5539, 5662, 5777, 5839, 5866, 5926, 5947, 6021, 6115, 6334, 6452, 6497, 6579, 6657, 6769, 6891, 6969, 7050, 7123, 7161, 7282, 7322, 7468, 7597, 7735, 7774, 7860, 7957, 8052, 8112, 8341, 8463, 8528, 8707, 8846, 8897, 8950, 9051, 9150, 9193, 9241, 9253, 9340, 9453, 9587, 9629, 9739, 9745 };

    long list_y1k_start[101] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 229, 234, 239, 244, 249, 254, 259, 264, 269, 274, 279, 284, 289, 294, 299, 304, 309, 314, 319, 324, 329, 334, 339, 344, 349, 354, 359, 364, 369, 374, 379, 384, 389, 394, 399, 404, 409, 414, 419, 424, 429, 434, 439, 444, 449, 454, 459, 464, 469, 474, 479, 484, 489, 494, 499 };

    long list_y1k_end[101] = { 4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79, 84, 89, 94, 99, 104, 109, 114, 119, 124, 129, 134, 139, 144, 149, 154, 159, 164, 169, 174, 179, 184, 189, 194, 199, 204, 209, 214, 219, 224, 228, 233, 238, 243, 248, 253, 258, 263, 268, 273, 278, 283, 288, 293, 298, 303, 308, 313, 318, 323, 328, 333, 338, 343, 348, 353, 358, 363, 368, 373, 378, 383, 388, 393, 398, 403, 408, 413, 418, 423, 428, 433, 438, 443, 448, 453, 458, 463, 468, 473, 478, 483, 488, 493, 498, 499 };

    long list_x10k_start[13] = { 0, 94, 185, 251, 354, 452, 498, 585, 666, 748, 834, 944, 1008 };

    long list_x10k_end[13] = { 93, 184, 252, 354, 453, 497, 585, 665, 748, 833, 944, 1007, 1044 };

    long list_y10k_start[13] = { 0, 42, 84, 126, 168, 210, 229, 271, 313, 355, 397, 439, 481 };

    long list_y10k_end[13] = { 41, 83, 125, 167, 209, 228, 270, 312, 354, 396, 438, 480, 499 };

    // the size of list_x_start_full would be the maximum SNP numbers
    long list_x_start_full[117904];


    long list_x10k_start_2[38] = { 0, 57, 83, 147, 153, 211, 265, 305, 312, 335, 364, 381, 397, 412, 423, 457, 467, 489, 498, 514, 530, 541, 558, 591, 598, 641, 651, 694, 725, 763, 786, 807, 878, 890, 921, 928, 966, 999 };
    long list_x10k_end_2[38] = { 55, 82, 146, 152, 210, 264, 304, 311, 334, 363, 380, 396, 411, 422, 456, 466, 488, 497, 513, 529, 540, 557, 590, 597, 640, 650, 693, 724, 762, 785, 806, 877, 889, 920, 927, 965, 998, 1044 };
    long list_y10k_start_2[38] = { 0, 28, 39, 70, 73, 101, 131, 151, 152, 161, 172, 175, 178, 183, 185, 212, 214, 228, 229, 237, 243, 248, 254, 273, 274, 307, 310, 327, 341, 367, 378, 385, 418, 419, 433, 434, 454, 479 };
    long list_y10k_end_2[38] = { 27, 38, 69, 72, 100, 130, 150, 151, 160, 171, 174, 177, 182, 184, 211, 213, 227, 228, 236, 242, 247, 253, 272, 273, 306, 309, 326, 340, 366, 377, 384, 417, 418, 432, 433, 453, 478, 499 };

    long list_x10k_start_3[500] = {0, 0, 0, 0, 0, 2, 2, 4, 6, 14, 14, 16, 16, 16, 18, 18, 18, 18, 22, 24, 24, 26, 26, 30, 30, 30, 30, 30, 54, 56, 58, 58, 58, 58, 58, 58, 62, 62, 66, 74, 74, 80, 82, 84, 86, 86, 86, 92, 92, 94, 96, 98, 104, 106, 108, 110, 110, 112, 114, 118, 120, 124, 124, 124, 126, 128, 128, 128, 130, 132, 138, 140, 140, 148, 148, 154, 154, 156, 158, 158, 160, 168, 170, 170, 172, 174, 174, 176, 178, 178, 180, 180, 182, 182, 182, 186, 186, 186, 190, 192, 192, 200, 200, 202, 206, 212, 212, 212, 214, 214, 216, 216, 220, 222, 224, 226, 226, 226, 226, 230, 232, 238, 238, 238, 238, 240, 240, 242, 244, 244, 244, 258, 258, 260, 264, 264, 268, 270, 272, 278, 280, 280, 280, 284, 284, 284, 286, 288, 290, 290, 290, 294, 302, 302, 306, 306, 308, 312, 314, 316, 318, 326, 332, 332, 336, 338, 338, 342, 342, 342, 342, 344, 356, 358, 360, 372, 374, 380, 388, 390, 392, 394, 396, 404, 406, 414, 414, 414, 416, 416, 418, 420, 420, 422, 422, 424, 424, 426, 428, 430, 430, 432, 434, 436, 436, 436, 436, 442, 442, 442, 442, 442, 448, 448, 458, 460, 462, 462, 462, 464, 464, 468, 468, 468, 470, 470, 470, 472, 472, 498, 498, 498, 498, 498, 498, 498, 498, 502, 504, 506, 510, 512, 514, 522, 524, 526, 528, 528, 532, 534, 534, 536, 536, 538, 550, 552, 554, 556, 556, 558, 562, 562, 564, 564, 564, 566, 568, 570, 570, 572, 572, 572, 574, 580, 590, 590, 592, 594, 594, 598, 600, 600, 600, 602, 602, 604, 604, 606, 606, 606, 606, 610, 610, 612, 612, 612, 612, 614, 614, 616, 616, 618, 618, 622, 622, 622, 624, 632, 632, 634, 644, 646, 652, 656, 658, 660, 664, 666, 670, 670, 670, 670, 672, 676, 676, 678, 678, 686, 686, 686, 686, 688, 690, 692, 694, 698, 698, 698, 704, 704, 706, 714, 714, 716, 718, 720, 722, 724, 726, 726, 730, 732, 732, 732, 734, 734, 736, 740, 742, 742, 742, 744, 744, 746, 746, 746, 748, 756, 758, 758, 758, 760, 760, 760, 762, 764, 764, 768, 776, 778, 780, 782, 786, 790, 790, 800, 804, 806, 806, 806, 808, 808, 808, 810, 812, 814, 818, 820, 820, 822, 828, 830, 830, 832, 838, 840, 840, 840, 844, 844, 846, 846, 848, 848, 850, 850, 854, 854, 874, 882, 884, 888, 890, 890, 892, 896, 896, 896, 896, 898, 902, 902, 904, 912, 920, 922, 924, 926, 928, 930, 930, 930, 932, 934, 940, 940, 942, 942, 942, 944, 946, 946, 946, 950, 956, 956, 960, 962, 962, 966, 966, 968, 968, 968, 970, 970, 970, 974, 974, 974, 974, 976, 976, 976, 976, 976, 980, 980, 984, 990, 992, 998, 998, 1000, 1000, 1010, 1010, 1012, 1014, 1016, 1016, 1018, 1018, 1018, 1018, 1018, 1018, 1018, 1018, 1018};
    long list_x1k_start_1[500] = {0, 70, 74, 84, 88, 116, 122, 130, 140, 212, 222, 226, 230, 242, 246, 250, 266, 272, 302, 308, 316, 324, 328, 352, 370, 386, 414, 428, 504, 522, 540, 560, 584, 586, 598, 618, 644, 646, 688, 748, 750, 804, 840, 852, 858, 858, 868, 916, 936, 950, 972, 990, 1038, 1056, 1086, 1096, 1104, 1118, 1140, 1180, 1204, 1234, 1244, 1250, 1254, 1270, 1272, 1284, 1298, 1312, 1366, 1384, 1396, 1476, 1478, 1524, 1532, 1550, 1554, 1558, 1580, 1652, 1664, 1670, 1698, 1702, 1710, 1720, 1746, 1750, 1754, 1760, 1770, 1778, 1780, 1818, 1820, 1822, 1868, 1882, 1890, 1966, 1974, 1986, 2034, 2090, 2098, 2104, 2120, 2120, 2146, 2148, 2186, 2212, 2216, 2242, 2242, 2250, 2252, 2280, 2300, 2350, 2352, 2356, 2356, 2380, 2386, 2400, 2430, 2432, 2436, 2548, 2556, 2564, 2616, 2616, 2640, 2672, 2678, 2734, 2754, 2756, 2758, 2788, 2794, 2802, 2810, 2854, 2856, 2858, 2864, 2922, 2988, 2994, 3016, 3028, 3036, 3088, 3094, 3118, 3158, 3220, 3260, 3282, 3312, 3330, 3332, 3352, 3360, 3364, 3364, 3390, 3482, 3506, 3522, 3640, 3656, 3716, 3788, 3804, 3820, 3840, 3852, 3954, 3956, 4020, 4022, 4022, 4022, 4026, 4030, 4054, 4060, 4068, 4080, 4092, 4100, 4128, 4160, 4164, 4172, 4192, 4200, 4216, 4220, 4228, 4228, 4266, 4268, 4272, 4274, 4274, 4322, 4328, 4394, 4424, 4438, 4444, 4446, 4458, 4486, 4514, 4516, 4524, 4532, 4538, 4538, 4544, 4610, 4718, 4718, 4726, 4730, 4744, 4774, 4788, 4798, 4850, 4854, 4874, 4914, 4948, 4956, 5036, 5052, 5068, 5086, 5088, 5138, 5146, 5158, 5170, 5184, 5196, 5304, 5328, 5336, 5354, 5356, 5370, 5402, 5408, 5420, 5428, 5436, 5462, 5476, 5492, 5500, 5508, 5516, 5518, 5526, 5590, 5670, 5670, 5694, 5716, 5734, 5756, 5766, 5776, 5776, 5802, 5804, 5808, 5808, 5832, 5832, 5838, 5844, 5874, 5874, 5890, 5892, 5894, 5896, 5898, 5904, 5914, 5916, 5938, 5948, 5980, 5990, 5990, 5994, 6068, 6070, 6088, 6170, 6196, 6250, 6288, 6308, 6316, 6364, 6394, 6412, 6418, 6418, 6420, 6430, 6462, 6466, 6478, 6482, 6534, 6544, 6548, 6558, 6566, 6594, 6612, 6632, 6662, 6662, 6668, 6728, 6740, 6752, 6832, 6832, 6842, 6870, 6872, 6900, 6924, 6930, 6942, 6980, 6996, 6998, 7010, 7016, 7024, 7050, 7080, 7088, 7094, 7100, 7112, 7118, 7122, 7128, 7136, 7150, 7236, 7250, 7250, 7256, 7266, 7274, 7276, 7298, 7326, 7328, 7356, 7428, 7438, 7462, 7488, 7508, 7556, 7566, 7654, 7688, 7698, 7698, 7700, 7714, 7716, 7728, 7732, 7746, 7770, 7800, 7818, 7818, 7836, 7868, 7892, 7896, 7902, 7946, 7968, 7974, 7980, 8016, 8022, 8036, 8052, 8056, 8066, 8088, 8092, 8126, 8130, 8264, 8346, 8348, 8402, 8414, 8426, 8432, 8472, 8478, 8484, 8484, 8504, 8546, 8546, 8550, 8636, 8706, 8726, 8754, 8772, 8806, 8818, 8820, 8820, 8824, 8840, 8886, 8890, 8908, 8912, 8912, 8920, 8942, 8952, 8954, 8990, 9046, 9054, 9080, 9094, 9098, 9134, 9136, 9148, 9152, 9154, 9162, 9166, 9172, 9200, 9206, 9208, 9208, 9214, 9214, 9214, 9218, 9220, 9254, 9258, 9280, 9332, 9352, 9392, 9398, 9416, 9420, 9502, 9508, 9518, 9550, 9556, 9564, 9580, 9586, 9590, 9604, 9638, 9664, 9676, 9676, 9676};

    long datatype;
    long data_version;
    long num_tag_snp;
    long num_sample;
    long num_target_snp;
    long window;
    long num_target_snp_k;

    long* list_x_start;
    long* list_x_end;
    long* list_y_start;
    long* list_y_end;

    long* list_nSNP_target;
    long* list_nSNP_tag;

    // SNP_Parameters(long _datatype, long = 0);

    SNP_Parameters(long _datatype, long _data_version, long _window, long num_target_snp_k);

    void setLists_newdata();
};

#endif