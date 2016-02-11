bool good_run(Int_t runNum, Int_t lumiNum) {
if(
        (runNum == 260532 && ((lumiNum >= 3 && lumiNum <= 8) || (lumiNum >= 10 && lumiNum <= 456) || (lumiNum >= 458 && lumiNum <= 746))) ||
        (runNum == 260533 && ((lumiNum >= 1 && lumiNum <= 14))) ||
        (runNum == 258713 && ((lumiNum >= 1 && lumiNum <= 161))) ||
        (runNum == 258712 && ((lumiNum >= 1 && lumiNum <= 524))) ||
        (runNum == 260536 && ((lumiNum >= 9 && lumiNum <= 37) || (lumiNum >= 45 && lumiNum <= 60) || (lumiNum >= 62 && lumiNum <= 193))) ||
        (runNum == 258714 && ((lumiNum >= 1 && lumiNum <= 67))) ||
        (runNum == 260534 && ((lumiNum >= 1 && lumiNum <= 375))) ||
        (runNum == 260431 && ((lumiNum >= 1 && lumiNum <= 447))) ||
        (runNum == 259637 && ((lumiNum >= 1 && lumiNum <= 72) || (lumiNum >= 75 && lumiNum <= 221))) ||
        (runNum == 258403 && ((lumiNum >= 1 && lumiNum <= 251))) ||
        (runNum == 258655 && ((lumiNum >= 60 && lumiNum <= 68))) ||
        (runNum == 258656 && ((lumiNum >= 1 && lumiNum <= 334))) ||
        (runNum == 258159 && ((lumiNum >= 1 && lumiNum <= 501))) ||
        (runNum == 258158 && ((lumiNum >= 1 && lumiNum <= 1088) || (lumiNum >= 1091 && lumiNum <= 1786))) ||
        (runNum == 258157 && ((lumiNum >= 1 && lumiNum <= 56))) ||
        (runNum == 256869 && ((lumiNum >= 1 && lumiNum <= 34))) ||
        (runNum == 256868 && ((lumiNum >= 5 && lumiNum <= 33) || (lumiNum >= 35 && lumiNum <= 200) || (lumiNum >= 202 && lumiNum <= 492))) ||
        (runNum == 256867 && ((lumiNum >= 1 && lumiNum <= 16) || (lumiNum >= 19 && lumiNum <= 94))) ||
        (runNum == 256866 && ((lumiNum >= 34 && lumiNum <= 47))) ||
        (runNum == 260538 && ((lumiNum >= 1 && lumiNum <= 284))) ||
        (runNum == 258434 && ((lumiNum >= 1 && lumiNum <= 453))) ||
        (runNum == 259721 && ((lumiNum >= 73 && lumiNum <= 99) || (lumiNum >= 102 && lumiNum <= 408))) ||
        (runNum == 257461 && ((lumiNum >= 44 && lumiNum <= 95))) ||
        (runNum == 258750 && ((lumiNum >= 1 && lumiNum <= 197))) ||
        (runNum == 258432 && ((lumiNum >= 1 && lumiNum <= 4))) ||
        (runNum == 259681 && ((lumiNum >= 64 && lumiNum <= 98))) ||
        (runNum == 259683 && ((lumiNum >= 5 && lumiNum <= 19) || (lumiNum >= 22 && lumiNum <= 23) || (lumiNum >= 25 && lumiNum <= 94))) ||
        (runNum == 259685 && ((lumiNum >= 1 && lumiNum <= 209) || (lumiNum >= 213 && lumiNum <= 240) || (lumiNum >= 242 && lumiNum <= 290) || (lumiNum >= 292 && lumiNum <= 445) || (lumiNum >= 447 && lumiNum <= 538) || (lumiNum >= 540 && lumiNum <= 544) || (lumiNum >= 546 && lumiNum <= 630))) ||
        (runNum == 259686 && ((lumiNum >= 1 && lumiNum <= 43) || (lumiNum >= 45 && lumiNum <= 47) || (lumiNum >= 49 && lumiNum <= 100) || (lumiNum >= 102 && lumiNum <= 108) || (lumiNum >= 110 && lumiNum <= 163) || (lumiNum >= 165 && lumiNum <= 245) || (lumiNum >= 248 && lumiNum <= 341))) ||
        (runNum == 254879 && ((lumiNum >= 52 && lumiNum <= 52) || (lumiNum >= 54 && lumiNum <= 140))) ||
        (runNum == 254231 && ((lumiNum >= 1 && lumiNum <= 24))) ||
        (runNum == 258213 && ((lumiNum >= 1 && lumiNum <= 165))) ||
        (runNum == 254232 && ((lumiNum >= 1 && lumiNum <= 81))) ||
        (runNum == 258694 && ((lumiNum >= 23 && lumiNum <= 199))) ||
        (runNum == 258129 && ((lumiNum >= 30 && lumiNum <= 124))) ||
        (runNum == 258427 && ((lumiNum >= 1 && lumiNum <= 107))) ||
        (runNum == 258426 && ((lumiNum >= 1 && lumiNum <= 10))) ||
        (runNum == 258425 && ((lumiNum >= 3 && lumiNum <= 136))) ||
        (runNum == 258428 && ((lumiNum >= 1 && lumiNum <= 159))) ||
        (runNum == 260593 && ((lumiNum >= 65 && lumiNum <= 401))) ||
        (runNum == 256843 && ((lumiNum >= 1 && lumiNum <= 204) || (lumiNum >= 207 && lumiNum <= 284) || (lumiNum >= 286 && lumiNum <= 378) || (lumiNum >= 380 && lumiNum <= 461) || (lumiNum >= 463 && lumiNum <= 587) || (lumiNum >= 598 && lumiNum <= 627) || (lumiNum >= 630 && lumiNum <= 661) || (lumiNum >= 1001 && lumiNum <= 1034) || (lumiNum >= 1036 && lumiNum <= 1081) || (lumiNum >= 1083 && lumiNum <= 1191) || (lumiNum >= 1193 && lumiNum <= 1193) || (lumiNum >= 1195 && lumiNum <= 1329) || (lumiNum >= 1331 && lumiNum <= 1332))) ||
        (runNum == 256842 && ((lumiNum >= 131 && lumiNum <= 132))) ||
        (runNum == 256926 && ((lumiNum >= 35 && lumiNum <= 50) || (lumiNum >= 53 && lumiNum <= 62) || (lumiNum >= 64 && lumiNum <= 65))) ||
        (runNum == 258136 && ((lumiNum >= 1 && lumiNum <= 60))) ||
        (runNum == 260424 && ((lumiNum >= 3 && lumiNum <= 12) || (lumiNum >= 15 && lumiNum <= 266) || (lumiNum >= 269 && lumiNum <= 672))) ||
        (runNum == 260425 && ((lumiNum >= 1 && lumiNum <= 18) || (lumiNum >= 21 && lumiNum <= 55) || (lumiNum >= 58 && lumiNum <= 256))) ||
        (runNum == 260426 && ((lumiNum >= 1 && lumiNum <= 52) || (lumiNum >= 55 && lumiNum <= 296) || (lumiNum >= 298 && lumiNum <= 307) || (lumiNum >= 310 && lumiNum <= 504))) ||
        (runNum == 257682 && ((lumiNum >= 66 && lumiNum <= 366))) ||
        (runNum == 258211 && ((lumiNum >= 43 && lumiNum <= 129))) ||
        (runNum == 256674 && ((lumiNum >= 1 && lumiNum <= 2))) ||
        (runNum == 256675 && ((lumiNum >= 1 && lumiNum <= 106) || (lumiNum >= 111 && lumiNum <= 164))) ||
        (runNum == 256676 && ((lumiNum >= 1 && lumiNum <= 160) || (lumiNum >= 162 && lumiNum <= 208))) ||
        (runNum == 256677 && ((lumiNum >= 1 && lumiNum <= 291) || (lumiNum >= 293 && lumiNum <= 390) || (lumiNum >= 392 && lumiNum <= 397) || (lumiNum >= 400 && lumiNum <= 455) || (lumiNum >= 457 && lumiNum <= 482))) ||
        (runNum == 256673 && ((lumiNum >= 55 && lumiNum <= 56))) ||
        (runNum == 258745 && ((lumiNum >= 1 && lumiNum <= 260))) ||
        (runNum == 258742 && ((lumiNum >= 2 && lumiNum <= 693))) ||
        (runNum == 258741 && ((lumiNum >= 22 && lumiNum <= 72))) ||
        (runNum == 257804 && ((lumiNum >= 1 && lumiNum <= 17))) ||
        (runNum == 257805 && ((lumiNum >= 1 && lumiNum <= 249))) ||
        (runNum == 258749 && ((lumiNum >= 1 && lumiNum <= 204) || (lumiNum >= 220 && lumiNum <= 604))) ||
        (runNum == 260373 && ((lumiNum >= 47 && lumiNum <= 370) || (lumiNum >= 373 && lumiNum <= 408))) ||
        (runNum == 258445 && ((lumiNum >= 1 && lumiNum <= 302))) ||
        (runNum == 258444 && ((lumiNum >= 1 && lumiNum <= 37))) ||
        (runNum == 258446 && ((lumiNum >= 1 && lumiNum <= 142))) ||
        (runNum == 258440 && ((lumiNum >= 1 && lumiNum <= 442) || (lumiNum >= 444 && lumiNum <= 732))) ||
        (runNum == 260576 && ((lumiNum >= 2 && lumiNum <= 88) || (lumiNum >= 90 && lumiNum <= 150))) ||
        (runNum == 260577 && ((lumiNum >= 1 && lumiNum <= 76))) ||
        (runNum == 260575 && ((lumiNum >= 1 && lumiNum <= 24))) ||
        (runNum == 258448 && ((lumiNum >= 2 && lumiNum <= 100) || (lumiNum >= 102 && lumiNum <= 731))) ||
        (runNum == 257751 && ((lumiNum >= 1 && lumiNum <= 463))) ||
        (runNum == 257599 && ((lumiNum >= 42 && lumiNum <= 118))) ||
        (runNum == 258214 && ((lumiNum >= 1 && lumiNum <= 217))) ||
        (runNum == 258215 && ((lumiNum >= 1 && lumiNum <= 6))) ||
        (runNum == 257613 && ((lumiNum >= 14 && lumiNum <= 1307))) ||
        (runNum == 257614 && ((lumiNum >= 1 && lumiNum <= 16))) ||
        (runNum == 257816 && ((lumiNum >= 1 && lumiNum <= 385))) ||
        (runNum == 257819 && ((lumiNum >= 1 && lumiNum <= 248))) ||
        (runNum == 259861 && ((lumiNum >= 1 && lumiNum <= 34) || (lumiNum >= 36 && lumiNum <= 38) || (lumiNum >= 40 && lumiNum <= 66) || (lumiNum >= 69 && lumiNum <= 77))) ||
        (runNum == 259862 && ((lumiNum >= 1 && lumiNum <= 13) || (lumiNum >= 16 && lumiNum <= 532))) ||
        (runNum == 256941 && ((lumiNum >= 1 && lumiNum <= 17) || (lumiNum >= 19 && lumiNum <= 29) || (lumiNum >= 103 && lumiNum <= 105) || (lumiNum >= 107 && lumiNum <= 126) || (lumiNum >= 129 && lumiNum <= 129) || (lumiNum >= 131 && lumiNum <= 168) || (lumiNum >= 170 && lumiNum <= 170) || (lumiNum >= 175 && lumiNum <= 290) || (lumiNum >= 293 && lumiNum <= 294))) ||
        (runNum == 258287 && ((lumiNum >= 45 && lumiNum <= 144) || (lumiNum >= 148 && lumiNum <= 227))) ||
        (runNum == 259809 && ((lumiNum >= 53 && lumiNum <= 222))) ||
        (runNum == 259810 && ((lumiNum >= 1 && lumiNum <= 113) || (lumiNum >= 116 && lumiNum <= 116))) ||
        (runNum == 259884 && ((lumiNum >= 73 && lumiNum <= 143) || (lumiNum >= 147 && lumiNum <= 155))) ||
        (runNum == 257968 && ((lumiNum >= 69 && lumiNum <= 326))) ||
        (runNum == 259818 && ((lumiNum >= 1 && lumiNum <= 160))) ||
        (runNum == 259817 && ((lumiNum >= 1 && lumiNum <= 5))) ||
        (runNum == 259811 && ((lumiNum >= 1 && lumiNum <= 47) || (lumiNum >= 50 && lumiNum <= 91))) ||
        (runNum == 257531 && ((lumiNum >= 5 && lumiNum <= 45) || (lumiNum >= 50 && lumiNum <= 143))) ||
        (runNum == 259813 && ((lumiNum >= 1 && lumiNum <= 10))) ||
        (runNum == 260627 && ((lumiNum >= 97 && lumiNum <= 611) || (lumiNum >= 613 && lumiNum <= 757) || (lumiNum >= 760 && lumiNum <= 788) || (lumiNum >= 791 && lumiNum <= 1051) || (lumiNum >= 1054 && lumiNum <= 1530) || (lumiNum >= 1533 && lumiNum <= 1845))) ||
        (runNum == 259891 && ((lumiNum >= 1 && lumiNum <= 108))) ||
        (runNum == 259890 && ((lumiNum >= 1 && lumiNum <= 34) || (lumiNum >= 37 && lumiNum <= 109))) ||
        (runNum == 254906 && ((lumiNum >= 1 && lumiNum <= 75))) ||
        (runNum == 254907 && ((lumiNum >= 1 && lumiNum <= 52))) ||
        (runNum == 258177 && ((lumiNum >= 1 && lumiNum <= 342) || (lumiNum >= 347 && lumiNum <= 724) || (lumiNum >= 755 && lumiNum <= 1939))) ||
        (runNum == 256801 && ((lumiNum >= 73 && lumiNum <= 263))) ||
        (runNum == 257735 && ((lumiNum >= 1 && lumiNum <= 15))) ||
        (runNum == 258702 && ((lumiNum >= 52 && lumiNum <= 402))) ||
        (runNum == 258703 && ((lumiNum >= 1 && lumiNum <= 389))) ||
        (runNum == 254914 && ((lumiNum >= 32 && lumiNum <= 32) || (lumiNum >= 34 && lumiNum <= 78))) ||
        (runNum == 258706 && ((lumiNum >= 1 && lumiNum <= 733))) ||
        (runNum == 260541 && ((lumiNum >= 1 && lumiNum <= 24))) ||
        (runNum == 258705 && ((lumiNum >= 1 && lumiNum <= 100))) ||
        (runNum == 259820 && ((lumiNum >= 1 && lumiNum <= 32) || (lumiNum >= 36 && lumiNum <= 161))) ||
        (runNum == 259821 && ((lumiNum >= 1 && lumiNum <= 75) || (lumiNum >= 78 && lumiNum <= 212))) ||
        (runNum == 259822 && ((lumiNum >= 1 && lumiNum <= 14) || (lumiNum >= 17 && lumiNum <= 464))) ||
        (runNum == 260427 && ((lumiNum >= 1 && lumiNum <= 198))) ||
        (runNum == 259626 && ((lumiNum >= 83 && lumiNum <= 106) || (lumiNum >= 108 && lumiNum <= 111) || (lumiNum >= 115 && lumiNum <= 166) || (lumiNum >= 169 && lumiNum <= 215) || (lumiNum >= 218 && lumiNum <= 437))) ||
        (runNum == 254852 && ((lumiNum >= 47 && lumiNum <= 94))) ||
        (runNum == 256630 && ((lumiNum >= 5 && lumiNum <= 26))) ||
        (runNum == 257969 && ((lumiNum >= 1 && lumiNum <= 634))) ||
        (runNum == 257645 && ((lumiNum >= 37 && lumiNum <= 73) || (lumiNum >= 75 && lumiNum <= 1096))) ||
        (runNum == 257723 && ((lumiNum >= 1 && lumiNum <= 1) || (lumiNum >= 3 && lumiNum <= 108) || (lumiNum >= 114 && lumiNum <= 148))) ||
        (runNum == 257722 && ((lumiNum >= 1 && lumiNum <= 19))) ||
        (runNum == 254790 && ((lumiNum >= 90 && lumiNum <= 90) || (lumiNum >= 93 && lumiNum <= 630) || (lumiNum >= 633 && lumiNum <= 697) || (lumiNum >= 701 && lumiNum <= 715) || (lumiNum >= 719 && lumiNum <= 784)))
  ) return true;
  else return false;

}
