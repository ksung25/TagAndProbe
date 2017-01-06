bool good_run(Int_t runNum, Int_t lumiNum) {
  if(
    (runNum == 273158 && ((lumiNum>=1 && lumiNum<=1279))) ||
    (runNum == 273302 && ((lumiNum>=1 && lumiNum<=459))) ||
    (runNum == 273402 && ((lumiNum>=100 && lumiNum<=292))) ||
    (runNum == 273403 && ((lumiNum>=1 && lumiNum<=53))) ||
    (runNum == 273404 && ((lumiNum>=1 && lumiNum<=18))) ||
    (runNum == 273405 && ((lumiNum>=2 && lumiNum<=25))) ||
    (runNum == 273406 && ((lumiNum>=1 && lumiNum<=112))) ||
    (runNum == 273408 && ((lumiNum>=1 && lumiNum<=6))) ||
    (runNum == 273409 && ((lumiNum>=1 && lumiNum<=309))) ||
    (runNum == 273410 && ((lumiNum>=1 && lumiNum<=90))) ||
    (runNum == 273411 && ((lumiNum>=1 && lumiNum<=29))) ||
    (runNum == 273425 && ((lumiNum>=62 && lumiNum<=352) || (lumiNum>=354 && lumiNum<=733))) ||
    (runNum == 273446 && ((lumiNum>=1 && lumiNum<=33))) ||
    (runNum == 273447 && ((lumiNum>=1 && lumiNum<=113) || (lumiNum>=115 && lumiNum<=412))) ||
    (runNum == 273448 && ((lumiNum>=1 && lumiNum<=391))) ||
    (runNum == 273449 && ((lumiNum>=1 && lumiNum<=214))) ||
    (runNum == 273450 && ((lumiNum>=1 && lumiNum<=214) || (lumiNum>=219 && lumiNum<=647))) ||
    (runNum == 273492 && ((lumiNum>=71 && lumiNum<=71) || (lumiNum>=73 && lumiNum<=282) || (lumiNum>=284 && lumiNum<=325) || (lumiNum>=327 && lumiNum<=338))) ||
    (runNum == 273493 && ((lumiNum>=1 && lumiNum<=233))) ||
    (runNum == 273494 && ((lumiNum>=1 && lumiNum<=192))) ||
    (runNum == 273502 && ((lumiNum>=73 && lumiNum<=256) || (lumiNum>=258 && lumiNum<=318) || (lumiNum>=320 && lumiNum<=813) || (lumiNum>=815 && lumiNum<=1064))) ||
    (runNum == 273503 && ((lumiNum>=1 && lumiNum<=598))) ||
    (runNum == 273554 && ((lumiNum>=77 && lumiNum<=437))) ||
    (runNum == 273555 && ((lumiNum>=1 && lumiNum<=173))) ||
    (runNum == 273725 && ((lumiNum>=83 && lumiNum<=252) || (lumiNum>=254 && lumiNum<=2545))) ||
    (runNum == 273728 && ((lumiNum>=1 && lumiNum<=100))) ||
    (runNum == 273730 && ((lumiNum>=1 && lumiNum<=1814) || (lumiNum>=1820 && lumiNum<=2126))) ||
    (runNum == 274094 && ((lumiNum>=108 && lumiNum<=332))) ||
    (runNum == 274146 && ((lumiNum>=1 && lumiNum<=67))) ||
    (runNum == 274157 && ((lumiNum>=105 && lumiNum<=534))) ||
    (runNum == 274159 && ((lumiNum>=1 && lumiNum<=43))) ||
    (runNum == 274160 && ((lumiNum>=1 && lumiNum<=207))) ||
    (runNum == 274161 && ((lumiNum>=1 && lumiNum<=516))) ||
    (runNum == 274172 && ((lumiNum>=31 && lumiNum<=95))) ||
    (runNum == 274198 && ((lumiNum>=81 && lumiNum<=191))) ||
    (runNum == 274199 && ((lumiNum>=1 && lumiNum<=623))) ||
    (runNum == 274200 && ((lumiNum>=1 && lumiNum<=678))) ||
    (runNum == 274240 && ((lumiNum>=1 && lumiNum<=40) || (lumiNum>=42 && lumiNum<=82))) ||
    (runNum == 274241 && ((lumiNum>=1 && lumiNum<=1152) || (lumiNum>=1161 && lumiNum<=1176))) ||
    (runNum == 274244 && ((lumiNum>=1 && lumiNum<=607))) ||
    (runNum == 274250 && ((lumiNum>=1 && lumiNum<=701))) ||
    (runNum == 274251 && ((lumiNum>=1 && lumiNum<=546))) ||
    (runNum == 274283 && ((lumiNum>=2 && lumiNum<=19))) ||
    (runNum == 274284 && ((lumiNum>=1 && lumiNum<=210))) ||
    (runNum == 274286 && ((lumiNum>=1 && lumiNum<=154))) ||
    (runNum == 274314 && ((lumiNum>=97 && lumiNum<=97) || (lumiNum>=99 && lumiNum<=158))) ||
    (runNum == 274315 && ((lumiNum>=1 && lumiNum<=424))) ||
    (runNum == 274316 && ((lumiNum>=1 && lumiNum<=959))) ||
    (runNum == 274317 && ((lumiNum>=1 && lumiNum<=3))) ||
    (runNum == 274319 && ((lumiNum>=1 && lumiNum<=225))) ||
    (runNum == 274335 && ((lumiNum>=60 && lumiNum<=1003))) ||
    (runNum == 274336 && ((lumiNum>=1 && lumiNum<=14))) ||
    (runNum == 274337 && ((lumiNum>=3 && lumiNum<=17))) ||
    (runNum == 274338 && ((lumiNum>=1 && lumiNum<=698))) ||
    (runNum == 274339 && ((lumiNum>=1 && lumiNum<=29) || (lumiNum>=31 && lumiNum<=31) || (lumiNum>=33 && lumiNum<=33) || (lumiNum>=35 && lumiNum<=93))) ||
    (runNum == 274344 && ((lumiNum>=1 && lumiNum<=632))) ||
    (runNum == 274345 && ((lumiNum>=1 && lumiNum<=170))) ||
    (runNum == 274382 && ((lumiNum>=94 && lumiNum<=144))) ||
    (runNum == 274387 && ((lumiNum>=88 && lumiNum<=439))) ||
    (runNum == 274388 && ((lumiNum>=1 && lumiNum<=1820))) ||
    (runNum == 274420 && ((lumiNum>=94 && lumiNum<=268))) ||
    (runNum == 274421 && ((lumiNum>=1 && lumiNum<=342))) ||
    (runNum == 274422 && ((lumiNum>=1 && lumiNum<=2207))) ||
    (runNum == 274440 && ((lumiNum>=92 && lumiNum<=493))) ||
    (runNum == 274441 && ((lumiNum>=1 && lumiNum<=431))) ||
    (runNum == 274442 && ((lumiNum>=1 && lumiNum<=752))) ||
    (runNum == 274954 && ((lumiNum>=37 && lumiNum<=37) || (lumiNum>=39 && lumiNum<=57))) ||
    (runNum == 274955 && ((lumiNum>=1 && lumiNum<=91))) ||
    (runNum == 274968 && ((lumiNum>=1 && lumiNum<=1192))) ||
    (runNum == 274969 && ((lumiNum>=1 && lumiNum<=1003))) ||
    (runNum == 274970 && ((lumiNum>=1 && lumiNum<=47))) ||
    (runNum == 274971 && ((lumiNum>=1 && lumiNum<=905))) ||
    (runNum == 274998 && ((lumiNum>=64 && lumiNum<=782))) ||
    (runNum == 274999 && ((lumiNum>=1 && lumiNum<=1241))) ||
    (runNum == 275000 && ((lumiNum>=1 && lumiNum<=136))) ||
    (runNum == 275001 && ((lumiNum>=1 && lumiNum<=1781) || (lumiNum>=1786 && lumiNum<=2061))) ||
    (runNum == 275059 && ((lumiNum>=78 && lumiNum<=81) || (lumiNum>=105 && lumiNum<=137))) ||
    (runNum == 275066 && ((lumiNum>=1 && lumiNum<=96))) ||
    (runNum == 275067 && ((lumiNum>=1 && lumiNum<=392))) ||
    (runNum == 275068 && ((lumiNum>=1 && lumiNum<=915))) ||
    (runNum == 275073 && ((lumiNum>=1 && lumiNum<=517))) ||
    (runNum == 275074 && ((lumiNum>=1 && lumiNum<=442) || (lumiNum>=444 && lumiNum<=647))) ||
    (runNum == 275124 && ((lumiNum>=106 && lumiNum<=106) || (lumiNum>=108 && lumiNum<=431))) ||
    (runNum == 275125 && ((lumiNum>=1 && lumiNum<=989))) ||
    (runNum == 275282 && ((lumiNum>=91 && lumiNum<=180))) ||
    (runNum == 275283 && ((lumiNum>=1 && lumiNum<=132))) ||
    (runNum == 275284 && ((lumiNum>=1 && lumiNum<=74))) ||
    (runNum == 275290 && ((lumiNum>=96 && lumiNum<=143))) ||
    (runNum == 275291 && ((lumiNum>=1 && lumiNum<=347))) ||
    (runNum == 275292 && ((lumiNum>=1 && lumiNum<=121))) ||
    (runNum == 275293 && ((lumiNum>=1 && lumiNum<=142) || (lumiNum>=144 && lumiNum<=201))) ||
    (runNum == 275309 && ((lumiNum>=55 && lumiNum<=617))) ||
    (runNum == 275310 && ((lumiNum>=1 && lumiNum<=1929))) ||
    (runNum == 275311 && ((lumiNum>=1 && lumiNum<=1253))) ||
    (runNum == 275319 && ((lumiNum>=141 && lumiNum<=282))) ||
    (runNum == 275337 && ((lumiNum>=1 && lumiNum<=427))) ||
    (runNum == 275338 && ((lumiNum>=1 && lumiNum<=520))) ||
    (runNum == 275344 && ((lumiNum>=76 && lumiNum<=356))) ||
    (runNum == 275345 && ((lumiNum>=1 && lumiNum<=353))) ||
    (runNum == 275370 && ((lumiNum>=81 && lumiNum<=365))) ||
    (runNum == 275371 && ((lumiNum>=1 && lumiNum<=22) || (lumiNum>=28 && lumiNum<=569))) ||
    (runNum == 275375 && ((lumiNum>=127 && lumiNum<=1449))) ||
    (runNum == 275376 && ((lumiNum>=1 && lumiNum<=2667) || (lumiNum>=2669 && lumiNum<=3096))) ||
    (runNum == 275657 && ((lumiNum>=1 && lumiNum<=105))) ||
    (runNum == 275658 && ((lumiNum>=1 && lumiNum<=337))) ||
    (runNum == 275659 && ((lumiNum>=1 && lumiNum<=17))) ||
    (runNum == 275761 && ((lumiNum>=1 && lumiNum<=9))) ||
    (runNum == 275767 && ((lumiNum>=1 && lumiNum<=4))) ||
    (runNum == 275772 && ((lumiNum>=1 && lumiNum<=56))) ||
    (runNum == 275773 && ((lumiNum>=1 && lumiNum<=7))) ||
    (runNum == 275774 && ((lumiNum>=1 && lumiNum<=311) || (lumiNum>=315 && lumiNum<=315))) ||
    (runNum == 275776 && ((lumiNum>=1 && lumiNum<=140))) ||
    (runNum == 275777 && ((lumiNum>=1 && lumiNum<=300))) ||
    (runNum == 275778 && ((lumiNum>=1 && lumiNum<=305))) ||
    (runNum == 275782 && ((lumiNum>=1 && lumiNum<=131) || (lumiNum>=133 && lumiNum<=762))) ||
    (runNum == 275832 && ((lumiNum>=1 && lumiNum<=367))) ||
    (runNum == 275833 && ((lumiNum>=1 && lumiNum<=53) || (lumiNum>=56 && lumiNum<=115) || (lumiNum>=117 && lumiNum<=251))) ||
    (runNum == 275834 && ((lumiNum>=1 && lumiNum<=297))) ||
    (runNum == 275835 && ((lumiNum>=1 && lumiNum<=13))) ||
    (runNum == 275836 && ((lumiNum>=1 && lumiNum<=429) || (lumiNum>=431 && lumiNum<=1163) || (lumiNum>=1166 && lumiNum<=1170) || (lumiNum>=1184 && lumiNum<=1293))) ||
    (runNum == 275837 && ((lumiNum>=1 && lumiNum<=186) || (lumiNum>=198 && lumiNum<=726))) ||
    (runNum == 275847 && ((lumiNum>=1 && lumiNum<=2263))) ||
    (runNum == 275886 && ((lumiNum>=73 && lumiNum<=109))) ||
    (runNum == 275890 && ((lumiNum>=1 && lumiNum<=1393))) ||
    (runNum == 275911 && ((lumiNum>=62 && lumiNum<=298) || (lumiNum>=300 && lumiNum<=354) || (lumiNum>=356 && lumiNum<=440))) ||
    (runNum == 275912 && ((lumiNum>=1 && lumiNum<=289))) ||
    (runNum == 275913 && ((lumiNum>=1 && lumiNum<=475))) ||
    (runNum == 275918 && ((lumiNum>=1 && lumiNum<=318) || (lumiNum>=348 && lumiNum<=361))) ||
    (runNum == 275920 && ((lumiNum>=5 && lumiNum<=463))) ||
    (runNum == 275921 && ((lumiNum>=1 && lumiNum<=2) || (lumiNum>=4 && lumiNum<=5) || (lumiNum>=17 && lumiNum<=20))) ||
    (runNum == 275923 && ((lumiNum>=3 && lumiNum<=53) || (lumiNum>=63 && lumiNum<=64) || (lumiNum>=66 && lumiNum<=126))) ||
    (runNum == 275931 && ((lumiNum>=1 && lumiNum<=14) || (lumiNum>=19 && lumiNum<=89))) ||
    (runNum == 275963 && ((lumiNum>=82 && lumiNum<=139) || (lumiNum>=141 && lumiNum<=172))) ||
    (runNum == 276092 && ((lumiNum>=74 && lumiNum<=149))) ||
    (runNum == 276097 && ((lumiNum>=1 && lumiNum<=507))) ||
    (runNum == 276242 && ((lumiNum>=1 && lumiNum<=7) || (lumiNum>=18 && lumiNum<=61) || (lumiNum>=72 && lumiNum<=1664))) ||
    (runNum == 276243 && ((lumiNum>=1 && lumiNum<=15) || (lumiNum>=18 && lumiNum<=480) || (lumiNum>=482 && lumiNum<=611))) ||
    (runNum == 276244 && ((lumiNum>=3 && lumiNum<=1202))) ||
    (runNum == 276282 && ((lumiNum>=75 && lumiNum<=534) || (lumiNum>=537 && lumiNum<=1142))) ||
    (runNum == 276283 && ((lumiNum>=3 && lumiNum<=1087))) ||
    (runNum == 276315 && ((lumiNum>=40 && lumiNum<=175) || (lumiNum>=178 && lumiNum<=217))) ||
    (runNum == 276317 && ((lumiNum>=3 && lumiNum<=138))) ||
    (runNum == 276318 && ((lumiNum>=3 && lumiNum<=103) || (lumiNum>=106 && lumiNum<=570))) ||
    (runNum == 276355 && ((lumiNum>=1 && lumiNum<=33))) ||
    (runNum == 276361 && ((lumiNum>=1 && lumiNum<=161) || (lumiNum>=169 && lumiNum<=208) || (lumiNum>=210 && lumiNum<=800) || (lumiNum>=802 && lumiNum<=833))) ||
    (runNum == 276363 && ((lumiNum>=1 && lumiNum<=140) || (lumiNum>=142 && lumiNum<=238) || (lumiNum>=242 && lumiNum<=1482))) ||
    (runNum == 276384 && ((lumiNum>=2 && lumiNum<=1117))) ||
    (runNum == 276437 && ((lumiNum>=63 && lumiNum<=224) || (lumiNum>=227 && lumiNum<=1074) || (lumiNum>=1076 && lumiNum<=2190))) ||
    (runNum == 276454 && ((lumiNum>=1 && lumiNum<=527))) ||
    (runNum == 276458 && ((lumiNum>=1 && lumiNum<=341))) ||
    (runNum == 276495 && ((lumiNum>=87 && lumiNum<=268))) ||
    (runNum == 276501 && ((lumiNum>=4 && lumiNum<=221) || (lumiNum>=223 && lumiNum<=2547))) ||
    (runNum == 276502 && ((lumiNum>=2 && lumiNum<=741))) ||
    (runNum == 276525 && ((lumiNum>=88 && lumiNum<=469) || (lumiNum>=471 && lumiNum<=1606) || (lumiNum>=1626 && lumiNum<=2893))) ||
    (runNum == 276527 && ((lumiNum>=1 && lumiNum<=214))) ||
    (runNum == 276528 && ((lumiNum>=4 && lumiNum<=394))) ||
    (runNum == 276542 && ((lumiNum>=74 && lumiNum<=857))) ||
    (runNum == 276543 && ((lumiNum>=1 && lumiNum<=638) || (lumiNum>=643 && lumiNum<=952))) ||
    (runNum == 276544 && ((lumiNum>=2 && lumiNum<=161))) ||
    (runNum == 276545 && ((lumiNum>=2 && lumiNum<=110) || (lumiNum>=117 && lumiNum<=213))) ||
    (runNum == 276581 && ((lumiNum>=79 && lumiNum<=444))) ||
    (runNum == 276582 && ((lumiNum>=1 && lumiNum<=871))) ||
    (runNum == 276583 && ((lumiNum>=1 && lumiNum<=52))) ||
    (runNum == 276584 && ((lumiNum>=1 && lumiNum<=2))) ||
    (runNum == 276585 && ((lumiNum>=1 && lumiNum<=238) || (lumiNum>=241 && lumiNum<=242) || (lumiNum>=245 && lumiNum<=246))) ||
    (runNum == 276586 && ((lumiNum>=2 && lumiNum<=658) || (lumiNum>=680 && lumiNum<=773))) ||
    (runNum == 276587 && ((lumiNum>=1 && lumiNum<=1006))) ||
    (runNum == 276653 && ((lumiNum>=72 && lumiNum<=550))) ||
    (runNum == 276655 && ((lumiNum>=1 && lumiNum<=593) || (lumiNum>=595 && lumiNum<=1106))) ||
    (runNum == 276659 && ((lumiNum>=1 && lumiNum<=127) || (lumiNum>=129 && lumiNum<=252))) ||
    (runNum == 276775 && ((lumiNum>=96 && lumiNum<=1260))) ||
    (runNum == 276776 && ((lumiNum>=1 && lumiNum<=1823))) ||
    (runNum == 276794 && ((lumiNum>=1 && lumiNum<=885))) ||
    (runNum == 276807 && ((lumiNum>=66 && lumiNum<=220))) ||
    (runNum == 276808 && ((lumiNum>=1 && lumiNum<=875))) ||
    (runNum == 276810 && ((lumiNum>=1 && lumiNum<=287))) ||
    (runNum == 276811 && ((lumiNum>=1 && lumiNum<=1270) || (lumiNum>=1272 && lumiNum<=2563))) ||
    (runNum == 276831 && ((lumiNum>=64 && lumiNum<=755) || (lumiNum>=761 && lumiNum<=2702))) ||
    (runNum == 276834 && ((lumiNum>=1 && lumiNum<=720))) ||
    (runNum == 276870 && ((lumiNum>=78 && lumiNum<=1354) || (lumiNum>=1356 && lumiNum<=3108) || (lumiNum>=3111 && lumiNum<=3258) || (lumiNum>=3260 && lumiNum<=3484))) ||
    (runNum == 276935 && ((lumiNum>=79 && lumiNum<=184) || (lumiNum>=186 && lumiNum<=838) || (lumiNum>=842 && lumiNum<=906))) ||
    (runNum == 276940 && ((lumiNum>=70 && lumiNum<=213))) ||
    (runNum == 276946 && ((lumiNum>=1 && lumiNum<=27))) ||
    (runNum == 276947 && ((lumiNum>=1 && lumiNum<=89) || (lumiNum>=91 && lumiNum<=126) || (lumiNum>=135 && lumiNum<=141))) ||
    (runNum == 276948 && ((lumiNum>=1 && lumiNum<=474))) ||
    (runNum == 276950 && ((lumiNum>=1 && lumiNum<=2353))) ||
    (runNum == 277069 && ((lumiNum>=81 && lumiNum<=265) || (lumiNum>=267 && lumiNum<=390))) ||
    (runNum == 277070 && ((lumiNum>=1 && lumiNum<=309) || (lumiNum>=311 && lumiNum<=1059))) ||
    (runNum == 277071 && ((lumiNum>=1 && lumiNum<=82) || (lumiNum>=90 && lumiNum<=178))) ||
    (runNum == 277072 && ((lumiNum>=1 && lumiNum<=253) || (lumiNum>=256 && lumiNum<=466))) ||
    (runNum == 277073 && ((lumiNum>=1 && lumiNum<=90))) ||
    (runNum == 277076 && ((lumiNum>=1 && lumiNum<=3) || (lumiNum>=5 && lumiNum<=7) || (lumiNum>=9 && lumiNum<=35) || (lumiNum>=38 && lumiNum<=1037))) ||
    (runNum == 277087 && ((lumiNum>=204 && lumiNum<=1191))) ||
    (runNum == 277094 && ((lumiNum>=1 && lumiNum<=161) || (lumiNum>=164 && lumiNum<=584))) ||
    (runNum == 277096 && ((lumiNum>=1 && lumiNum<=1309) || (lumiNum>=1311 && lumiNum<=2086))) ||
    (runNum == 277112 && ((lumiNum>=1 && lumiNum<=155))) ||
    (runNum == 277126 && ((lumiNum>=42 && lumiNum<=59))) ||
    (runNum == 277127 && ((lumiNum>=1 && lumiNum<=438) || (lumiNum>=440 && lumiNum<=902))) ||
    (runNum == 277148 && ((lumiNum>=83 && lumiNum<=190) || (lumiNum>=193 && lumiNum<=700))) ||
    (runNum == 277166 && ((lumiNum>=77 && lumiNum<=186) || (lumiNum>=188 && lumiNum<=431))) ||
    (runNum == 277168 && ((lumiNum>=1 && lumiNum<=1708) || (lumiNum>=1711 && lumiNum<=1822) || (lumiNum>=1824 && lumiNum<=2223))) ||
    (runNum == 277180 && ((lumiNum>=88 && lumiNum<=228))) ||
    (runNum == 277194 && ((lumiNum>=113 && lumiNum<=139) || (lumiNum>=144 && lumiNum<=497) || (lumiNum>=500 && lumiNum<=1115) || (lumiNum>=1117 && lumiNum<=1312) || (lumiNum>=1320 && lumiNum<=1749) || (lumiNum>=1754 && lumiNum<=2067) || (lumiNum>=2070 && lumiNum<=2070))) ||
    (runNum == 277305 && ((lumiNum>=62 && lumiNum<=744))) ||
    (runNum == 277420 && ((lumiNum>=84 && lumiNum<=84) || (lumiNum>=86 && lumiNum<=291) || (lumiNum>=293 && lumiNum<=346))) ||
    (runNum == 277981 && ((lumiNum>=82 && lumiNum<=83) || (lumiNum>=85 && lumiNum<=163))) ||
    (runNum == 277991 && ((lumiNum>=1 && lumiNum<=98))) ||
    (runNum == 277992 && ((lumiNum>=1 && lumiNum<=260) || (lumiNum>=262 && lumiNum<=312))) ||
    (runNum == 278017 && ((lumiNum>=77 && lumiNum<=97) || (lumiNum>=99 && lumiNum<=213) || (lumiNum>=215 && lumiNum<=512) || (lumiNum>=514 && lumiNum<=589))) ||
    (runNum == 278018 && ((lumiNum>=1 && lumiNum<=263) || (lumiNum>=265 && lumiNum<=422) || (lumiNum>=424 && lumiNum<=615) || (lumiNum>=617 && lumiNum<=627) || (lumiNum>=642 && lumiNum<=1011) || (lumiNum>=1020 && lumiNum<=1181))) ||
    (runNum == 278167 && ((lumiNum>=87 && lumiNum<=394) || (lumiNum>=397 && lumiNum<=1153) || (lumiNum>=1155 && lumiNum<=1660) || (lumiNum>=1662 && lumiNum<=1707) || (lumiNum>=1709 && lumiNum<=2258))) ||
    (runNum == 278175 && ((lumiNum>=1 && lumiNum<=88))) ||
    (runNum == 278193 && ((lumiNum>=77 && lumiNum<=231))) ||
    (runNum == 278239 && ((lumiNum>=76 && lumiNum<=339) || (lumiNum>=341 && lumiNum<=558) || (lumiNum>=560 && lumiNum<=740))) ||
    (runNum == 278240 && ((lumiNum>=1 && lumiNum<=64) || (lumiNum>=70 && lumiNum<=113) || (lumiNum>=115 && lumiNum<=1121) || (lumiNum>=1123 && lumiNum<=1296) || (lumiNum>=1299 && lumiNum<=1309))) ||
    (runNum == 278273 && ((lumiNum>=75 && lumiNum<=110))) ||
    (runNum == 278274 && ((lumiNum>=1 && lumiNum<=18) || (lumiNum>=20 && lumiNum<=85))) ||
    (runNum == 278288 && ((lumiNum>=67 && lumiNum<=81))) ||
    (runNum == 278289 && ((lumiNum>=1 && lumiNum<=42) || (lumiNum>=44 && lumiNum<=52))) ||
    (runNum == 278290 && ((lumiNum>=1 && lumiNum<=11))) ||
    (runNum == 278308 && ((lumiNum>=87 && lumiNum<=216) || (lumiNum>=219 && lumiNum<=587) || (lumiNum>=589 && lumiNum<=680) || (lumiNum>=683 && lumiNum<=1200) || (lumiNum>=1217 && lumiNum<=1410) || (lumiNum>=1413 && lumiNum<=1848) || (lumiNum>=1880 && lumiNum<=1880))) ||
    (runNum == 278310 && ((lumiNum>=1 && lumiNum<=32) || (lumiNum>=34 && lumiNum<=709))) ||
    (runNum == 278315 && ((lumiNum>=73 && lumiNum<=254) || (lumiNum>=256 && lumiNum<=661) || (lumiNum>=663 && lumiNum<=767))) ||
    (runNum == 278345 && ((lumiNum>=84 && lumiNum<=500) || (lumiNum>=503 && lumiNum<=831))) ||
    (runNum == 278346 && ((lumiNum>=1 && lumiNum<=117))) ||
    (runNum == 278349 && ((lumiNum>=1 && lumiNum<=401) || (lumiNum>=403 && lumiNum<=612) || (lumiNum>=632 && lumiNum<=633))) ||
    (runNum == 278366 && ((lumiNum>=1 && lumiNum<=453))) ||
    (runNum == 278406 && ((lumiNum>=85 && lumiNum<=360) || (lumiNum>=362 && lumiNum<=1682))) ||
    (runNum == 278509 && ((lumiNum>=91 && lumiNum<=1557))) ||
    (runNum == 278769 && ((lumiNum>=75 && lumiNum<=104))) ||
    (runNum == 278770 && ((lumiNum>=1 && lumiNum<=767))) ||
    (runNum == 278801 && ((lumiNum>=48 && lumiNum<=85))) ||
    (runNum == 278802 && ((lumiNum>=1 && lumiNum<=17))) ||
    (runNum == 278803 && ((lumiNum>=1 && lumiNum<=87) || (lumiNum>=91 && lumiNum<=133) || (lumiNum>=135 && lumiNum<=297) || (lumiNum>=299 && lumiNum<=323))) ||
    (runNum == 278804 && ((lumiNum>=1 && lumiNum<=4))) ||
    (runNum == 278805 && ((lumiNum>=3 && lumiNum<=26) || (lumiNum>=30 && lumiNum<=167) || (lumiNum>=170 && lumiNum<=193) || (lumiNum>=196 && lumiNum<=280) || (lumiNum>=283 && lumiNum<=284) || (lumiNum>=288 && lumiNum<=288))) ||
    (runNum == 278808 && ((lumiNum>=1 && lumiNum<=445) || (lumiNum>=447 && lumiNum<=462) || (lumiNum>=464 && lumiNum<=1793))) ||
    (runNum == 278820 && ((lumiNum>=17 && lumiNum<=1533))) ||
    (runNum == 278822 && ((lumiNum>=1 && lumiNum<=1627))) ||
    (runNum == 278873 && ((lumiNum>=70 && lumiNum<=129))) ||
    (runNum == 278874 && ((lumiNum>=1 && lumiNum<=273) || (lumiNum>=275 && lumiNum<=478))) ||
    (runNum == 278875 && ((lumiNum>=1 && lumiNum<=210) || (lumiNum>=212 && lumiNum<=834))) ||
    (runNum == 278923 && ((lumiNum>=55 && lumiNum<=467))) ||
    (runNum == 278957 && ((lumiNum>=79 && lumiNum<=227))) ||
    (runNum == 278962 && ((lumiNum>=68 && lumiNum<=408))) ||
    (runNum == 278963 && ((lumiNum>=1 && lumiNum<=23) || (lumiNum>=25 && lumiNum<=175))) ||
    (runNum == 278969 && ((lumiNum>=70 && lumiNum<=511) || (lumiNum>=514 && lumiNum<=1051) || (lumiNum>=1053 && lumiNum<=1291) || (lumiNum>=1293 && lumiNum<=1397) || (lumiNum>=1399 && lumiNum<=1460))) ||
    (runNum == 278975 && ((lumiNum>=1 && lumiNum<=475) || (lumiNum>=477 && lumiNum<=745) || (lumiNum>=747 && lumiNum<=850))) ||
    (runNum == 278976 && ((lumiNum>=1 && lumiNum<=20))) ||
    (runNum == 278986 && ((lumiNum>=71 && lumiNum<=199))) ||
    (runNum == 279024 && ((lumiNum>=82 && lumiNum<=382))) ||
    (runNum == 279029 && ((lumiNum>=1 && lumiNum<=260))) ||
    (runNum == 279071 && ((lumiNum>=71 && lumiNum<=244))) ||
    (runNum == 279080 && ((lumiNum>=68 && lumiNum<=224))) ||
    (runNum == 279115 && ((lumiNum>=118 && lumiNum<=524))) ||
    (runNum == 279116 && ((lumiNum>=38 && lumiNum<=485))) ||
    (runNum == 279479 && ((lumiNum>=86 && lumiNum<=190))) ||
    (runNum == 279588 && ((lumiNum>=100 && lumiNum<=1259))) ||
    (runNum == 279653 && ((lumiNum>=77 && lumiNum<=77) || (lumiNum>=82 && lumiNum<=261))) ||
    (runNum == 279654 && ((lumiNum>=1 && lumiNum<=108) || (lumiNum>=110 && lumiNum<=1231) || (lumiNum>=1285 && lumiNum<=1299))) ||
    (runNum == 279656 && ((lumiNum>=1 && lumiNum<=43))) ||
    (runNum == 279658 && ((lumiNum>=1 && lumiNum<=689) || (lumiNum>=691 && lumiNum<=713))) ||
    (runNum == 279667 && ((lumiNum>=68 && lumiNum<=1033))) ||
    (runNum == 279681 && ((lumiNum>=77 && lumiNum<=104))) ||
    (runNum == 279682 && ((lumiNum>=1 && lumiNum<=29) || (lumiNum>=33 && lumiNum<=34) || (lumiNum>=37 && lumiNum<=38))) ||
    (runNum == 279683 && ((lumiNum>=1 && lumiNum<=26))) ||
    (runNum == 279684 && ((lumiNum>=1 && lumiNum<=22))) ||
    (runNum == 279685 && ((lumiNum>=1 && lumiNum<=93) || (lumiNum>=95 && lumiNum<=209))) ||
    (runNum == 279691 && ((lumiNum>=71 && lumiNum<=113))) ||
    (runNum == 279694 && ((lumiNum>=1 && lumiNum<=2235))) ||
    (runNum == 279715 && ((lumiNum>=71 && lumiNum<=474) || (lumiNum>=476 && lumiNum<=477) || (lumiNum>=480 && lumiNum<=480) || (lumiNum>=511 && lumiNum<=511) || (lumiNum>=523 && lumiNum<=691))) ||
    (runNum == 279716 && ((lumiNum>=1 && lumiNum<=860) || (lumiNum>=875 && lumiNum<=1528) || (lumiNum>=1530 && lumiNum<=1653))) ||
    (runNum == 279760 && ((lumiNum>=68 && lumiNum<=578) || (lumiNum>=585 && lumiNum<=728))) ||
    (runNum == 279766 && ((lumiNum>=1 && lumiNum<=1689))) ||
    (runNum == 279767 && ((lumiNum>=1 && lumiNum<=776))) ||
    (runNum == 279794 && ((lumiNum>=77 && lumiNum<=1100))) ||
    (runNum == 279823 && ((lumiNum>=61 && lumiNum<=395))) ||
    (runNum == 279841 && ((lumiNum>=75 && lumiNum<=398) || (lumiNum>=407 && lumiNum<=2122))) ||
    (runNum == 279844 && ((lumiNum>=72 && lumiNum<=295))) ||
    (runNum == 279887 && ((lumiNum>=79 && lumiNum<=221) || (lumiNum>=225 && lumiNum<=397))) ||
    (runNum == 279931 && ((lumiNum>=84 && lumiNum<=628) || (lumiNum>=630 && lumiNum<=743) || (lumiNum>=746 && lumiNum<=801) || (lumiNum>=803 && lumiNum<=1043) || (lumiNum>=1045 && lumiNum<=3022))) ||
    (runNum == 279966 && ((lumiNum>=79 && lumiNum<=441))) ||
    (runNum == 279975 && ((lumiNum>=70 && lumiNum<=190) || (lumiNum>=192 && lumiNum<=253) || (lumiNum>=256 && lumiNum<=281) || (lumiNum>=283 && lumiNum<=709) || (lumiNum>=734 && lumiNum<=1121))) ||
    (runNum == 279993 && ((lumiNum>=85 && lumiNum<=156))) ||
    (runNum == 279994 && ((lumiNum>=1 && lumiNum<=47))) ||
    (runNum == 280013 && ((lumiNum>=1 && lumiNum<=25))) ||
    (runNum == 280015 && ((lumiNum>=1 && lumiNum<=39) || (lumiNum>=41 && lumiNum<=56) || (lumiNum>=59 && lumiNum<=554) || (lumiNum>=560 && lumiNum<=580))) ||
    (runNum == 280016 && ((lumiNum>=1 && lumiNum<=149))) ||
    (runNum == 280017 && ((lumiNum>=1 && lumiNum<=608))) ||
    (runNum == 280018 && ((lumiNum>=1 && lumiNum<=1281))) ||
    (runNum == 280020 && ((lumiNum>=1 && lumiNum<=45))) ||
    (runNum == 280024 && ((lumiNum>=1 && lumiNum<=427))) ||
    (runNum == 280187 && ((lumiNum>=4 && lumiNum<=60))) ||
    (runNum == 280188 && ((lumiNum>=1 && lumiNum<=245))) ||
    (runNum == 280191 && ((lumiNum>=1 && lumiNum<=781) || (lumiNum>=783 && lumiNum<=866) || (lumiNum>=869 && lumiNum<=900))) ||
    (runNum == 280194 && ((lumiNum>=1 && lumiNum<=238))) ||
    (runNum == 280242 && ((lumiNum>=1 && lumiNum<=411) || (lumiNum>=414 && lumiNum<=627))) ||
    (runNum == 280249 && ((lumiNum>=1 && lumiNum<=486) || (lumiNum>=488 && lumiNum<=1433))) ||
    (runNum == 280251 && ((lumiNum>=1 && lumiNum<=165) || (lumiNum>=167 && lumiNum<=372))) ||
    (runNum == 280327 && ((lumiNum>=49 && lumiNum<=85))) ||
    (runNum == 280330 && ((lumiNum>=1 && lumiNum<=857))) ||
    (runNum == 280349 && ((lumiNum>=1 && lumiNum<=247) || (lumiNum>=252 && lumiNum<=623) || (lumiNum>=626 && lumiNum<=626))) ||
    (runNum == 280363 && ((lumiNum>=1 && lumiNum<=359))) ||
    (runNum == 280364 && ((lumiNum>=1 && lumiNum<=370) || (lumiNum>=372 && lumiNum<=617) || (lumiNum>=619 && lumiNum<=619) || (lumiNum>=621 && lumiNum<=1090) || (lumiNum>=1102 && lumiNum<=1363))) ||
    (runNum == 280383 && ((lumiNum>=64 && lumiNum<=65))) ||
    (runNum == 280384 && ((lumiNum>=2 && lumiNum<=34))) ||
    (runNum == 280385 && ((lumiNum>=1 && lumiNum<=519) || (lumiNum>=523 && lumiNum<=569) || (lumiNum>=574 && lumiNum<=1187) || (lumiNum>=1189 && lumiNum<=1533) || (lumiNum>=1536 && lumiNum<=2022))) ||
    (runNum == 281613 && ((lumiNum>=101 && lumiNum<=128) || (lumiNum>=130 && lumiNum<=130) || (lumiNum>=133 && lumiNum<=133) || (lumiNum>=135 && lumiNum<=139) || (lumiNum>=143 && lumiNum<=256) || (lumiNum>=258 && lumiNum<=903))) ||
    (runNum == 281639 && ((lumiNum>=1 && lumiNum<=132))) ||
    (runNum == 281641 && ((lumiNum>=1 && lumiNum<=319))) ||
    (runNum == 281693 && ((lumiNum>=1 && lumiNum<=2191))) ||
    (runNum == 281707 && ((lumiNum>=99 && lumiNum<=982) || (lumiNum>=1000 && lumiNum<=1065))) ||
    (runNum == 281726 && ((lumiNum>=1 && lumiNum<=288))) ||
    (runNum == 281727 && ((lumiNum>=1 && lumiNum<=1605))) ||
    (runNum == 281797 && ((lumiNum>=125 && lumiNum<=2176))) ||
    (runNum == 281975 && ((lumiNum>=1 && lumiNum<=215))) ||
    (runNum == 281976 && ((lumiNum>=1 && lumiNum<=2166))) ||
    (runNum == 282033 && ((lumiNum>=82 && lumiNum<=117))) ||
    (runNum == 282034 && ((lumiNum>=1 && lumiNum<=33))) ||
    (runNum == 282035 && ((lumiNum>=1 && lumiNum<=40))) ||
    (runNum == 282037 && ((lumiNum>=1 && lumiNum<=457) || (lumiNum>=459 && lumiNum<=1862))) ||
    (runNum == 282092 && ((lumiNum>=92 && lumiNum<=222) || (lumiNum>=624 && lumiNum<=2276))) ||
    (runNum == 282708 && ((lumiNum>=1 && lumiNum<=8))) ||
    (runNum == 282710 && ((lumiNum>=1 && lumiNum<=2) || (lumiNum>=8 && lumiNum<=8))) ||
    (runNum == 282712 && ((lumiNum>=1 && lumiNum<=1) || (lumiNum>=10 && lumiNum<=68))) ||
    (runNum == 282730 && ((lumiNum>=89 && lumiNum<=164))) ||
    (runNum == 282731 && ((lumiNum>=1 && lumiNum<=172))) ||
    (runNum == 282732 && ((lumiNum>=1 && lumiNum<=69))) ||
    (runNum == 282733 && ((lumiNum>=1 && lumiNum<=177))) ||
    (runNum == 282734 && ((lumiNum>=1 && lumiNum<=327))) ||
    (runNum == 282735 && ((lumiNum>=1 && lumiNum<=642) || (lumiNum>=645 && lumiNum<=1232) || (lumiNum>=1235 && lumiNum<=1823))) ||
    (runNum == 282800 && ((lumiNum>=1 && lumiNum<=377))) ||
    (runNum == 282807 && ((lumiNum>=1 && lumiNum<=326))) ||
    (runNum == 282814 && ((lumiNum>=1 && lumiNum<=1843))) ||
    (runNum == 282842 && ((lumiNum>=1 && lumiNum<=80))) ||
    (runNum == 282917 && ((lumiNum>=117 && lumiNum<=157) || (lumiNum>=159 && lumiNum<=191))) ||
    (runNum == 282918 && ((lumiNum>=1 && lumiNum<=51))) ||
    (runNum == 282919 && ((lumiNum>=1 && lumiNum<=243))) ||
    (runNum == 282922 && ((lumiNum>=1 && lumiNum<=131))) ||
    (runNum == 282923 && ((lumiNum>=1 && lumiNum<=17) || (lumiNum>=19 && lumiNum<=30) || (lumiNum>=32 && lumiNum<=36) || (lumiNum>=38 && lumiNum<=39) || (lumiNum>=41 && lumiNum<=86) || (lumiNum>=88 && lumiNum<=224))) ||
    (runNum == 283042 && ((lumiNum>=1 && lumiNum<=6))) ||
    (runNum == 283043 && ((lumiNum>=1 && lumiNum<=105) || (lumiNum>=108 && lumiNum<=519))) ||
    (runNum == 283049 && ((lumiNum>=82 && lumiNum<=93))) ||
    (runNum == 283050 && ((lumiNum>=1 && lumiNum<=212))) ||
    (runNum == 283052 && ((lumiNum>=1 && lumiNum<=111))) ||
    (runNum == 283059 && ((lumiNum>=1 && lumiNum<=125) || (lumiNum>=127 && lumiNum<=451))) ||
    (runNum == 283270 && ((lumiNum>=76 && lumiNum<=573) || (lumiNum>=576 && lumiNum<=1502) || (lumiNum>=1504 && lumiNum<=1888) || (lumiNum>=1890 && lumiNum<=1912))) ||
    (runNum == 283283 && ((lumiNum>=4 && lumiNum<=1668) || (lumiNum>=1670 && lumiNum<=1748))) ||
    (runNum == 283305 && ((lumiNum>=79 && lumiNum<=85))) ||
    (runNum == 283306 && ((lumiNum>=1 && lumiNum<=289))) ||
    (runNum == 283307 && ((lumiNum>=1 && lumiNum<=153) || (lumiNum>=156 && lumiNum<=456))) ||
    (runNum == 283308 && ((lumiNum>=1 && lumiNum<=547) || (lumiNum>=549 && lumiNum<=571) || (lumiNum>=573 && lumiNum<=895) || (lumiNum>=897 && lumiNum<=948))) ||
    (runNum == 283353 && ((lumiNum>=80 && lumiNum<=822))) ||
    (runNum == 283358 && ((lumiNum>=1 && lumiNum<=243) || (lumiNum>=245 && lumiNum<=981))) ||
    (runNum == 283359 && ((lumiNum>=1 && lumiNum<=428))) ||
    (runNum == 283407 && ((lumiNum>=82 && lumiNum<=114))) ||
    (runNum == 283408 && ((lumiNum>=1 && lumiNum<=27) || (lumiNum>=29 && lumiNum<=2088) || (lumiNum>=2098 && lumiNum<=2125) || (lumiNum>=2203 && lumiNum<=2416) || (lumiNum>=2528 && lumiNum<=2542))) ||
    (runNum == 283416 && ((lumiNum>=49 && lumiNum<=151) || (lumiNum>=154 && lumiNum<=245))) ||
    (runNum == 283453 && ((lumiNum>=83 && lumiNum<=537))) ||
    (runNum == 283469 && ((lumiNum>=74 && lumiNum<=74))) ||
    (runNum == 283478 && ((lumiNum>=76 && lumiNum<=303) || (lumiNum>=324 && lumiNum<=969))) ||
    (runNum == 283548 && ((lumiNum>=145 && lumiNum<=288))) ||
    (runNum == 283680 && ((lumiNum>=1 && lumiNum<=81))) ||
    (runNum == 283681 && ((lumiNum>=1 && lumiNum<=17))) ||
    (runNum == 283682 && ((lumiNum>=1 && lumiNum<=384))) ||
    (runNum == 283685 && ((lumiNum>=1 && lumiNum<=314))) ||
    (runNum == 283820 && ((lumiNum>=67 && lumiNum<=1548))) ||
    (runNum == 283830 && ((lumiNum>=1 && lumiNum<=722))) ||
    (runNum == 283834 && ((lumiNum>=1 && lumiNum<=67) || (lumiNum>=69 && lumiNum<=82))) ||
    (runNum == 283835 && ((lumiNum>=1 && lumiNum<=14) || (lumiNum>=16 && lumiNum<=112))) ||
    (runNum == 283865 && ((lumiNum>=1 && lumiNum<=1177))) ||
    (runNum == 283876 && ((lumiNum>=65 && lumiNum<=211) || (lumiNum>=215 && lumiNum<=724))) ||
    (runNum == 283877 && ((lumiNum>=1 && lumiNum<=1496))) ||
    (runNum == 283884 && ((lumiNum>=349 && lumiNum<=504) || (lumiNum>=509 && lumiNum<=756))) ||
    (runNum == 283885 && ((lumiNum>=1 && lumiNum<=1723))) ||
    (runNum == 283933 && ((lumiNum>=88 && lumiNum<=232))) ||
    (runNum == 283934 && ((lumiNum>=1 && lumiNum<=784) || (lumiNum>=793 && lumiNum<=870) || (lumiNum>=875 && lumiNum<=1245) || (lumiNum>=1267 && lumiNum<=1291))) ||
    (runNum == 283946 && ((lumiNum>=85 && lumiNum<=1448) || (lumiNum>=1450 && lumiNum<=1462))) ||
    (runNum == 283964 && ((lumiNum>=1 && lumiNum<=388))) ||
    (runNum == 284006 && ((lumiNum>=73 && lumiNum<=390))) ||
    (runNum == 284014 && ((lumiNum>=1 && lumiNum<=266))) ||
    (runNum == 284025 && ((lumiNum>=110 && lumiNum<=157))) ||
    (runNum == 284029 && ((lumiNum>=1 && lumiNum<=112))) ||
    (runNum == 284035 && ((lumiNum>=1 && lumiNum<=360))) ||
    (runNum == 284036 && ((lumiNum>=1 && lumiNum<=140) || (lumiNum>=143 && lumiNum<=348))) ||
    (runNum == 284037 && ((lumiNum>=1 && lumiNum<=340))) ||
    (runNum == 284038 && ((lumiNum>=1 && lumiNum<=55))) ||
    (runNum == 284039 && ((lumiNum>=1 && lumiNum<=30))) ||
    (runNum == 284040 && ((lumiNum>=1 && lumiNum<=33))) ||
    (runNum == 284041 && ((lumiNum>=1 && lumiNum<=44))) ||
    (runNum == 284042 && ((lumiNum>=1 && lumiNum<=129))) ||
    (runNum == 284043 && ((lumiNum>=1 && lumiNum<=205) || (lumiNum>=210 && lumiNum<=224))) ||
    (runNum == 284044 && ((lumiNum>=1 && lumiNum<=30)))
  ) return true;
  else return false;

}
