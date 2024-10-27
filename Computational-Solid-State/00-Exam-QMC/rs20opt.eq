 nproc =            1
 enter run id
   ndim 2                                                                        
   type up 13 4.0000000000000001E-002 rs20opt.up.x                               
   type down 13 4.0000000000000001E-002 rs20opt.down.x                           
   v2 up up rs20opt.v 1                                                          
 routinename  fitpn                                          
 # param            0
          64  punti in spazio k
   v2 down down rs20opt.v 1                                                      
   v2 up down rs20opt.v 1                                                        
   v0 -14.275916300392243                                                        
   u2 up up rs20opt.u                                                            
 routinename  fitpn                                          
 # param            0
          64  punti in spazio k
   u2 down down rs20opt.u                                                        
   u2 up down rs20opt.u                                                          
   pbc 9.0377767727099023 9.0377767727099023                                     
   kspace rs20opt.k                                                              
   u3 up up rs20opt.u3                                                           
 routinename  csi                                            
 # param            3
   u3 up down rs20opt.u3                                                         
   u3 down down rs20opt.u3                                                       
   rhok up                                                                       
   rhok down                                                                     
   backflow up up rs20opt.b                                                      
 routinename  csi                                            
 # param            3
   backflow up down rs20opt.b                                                    
   backflow down down rs20opt.b                                                  
 seed         3186        1495        2205        5111
 n_props_in_stack =          270
 dmc 40 150 0.1 20 -0.297                                                        
 ===>> dmc block            1
 -0.29886688499E+00 0.308E+04 -0.29886688499E+00 0.00E+00  elocal              
  0.97280495253E+00 0.308E+04  0.97280495253E+00 0.00E+00  acc.rate            
 -0.37557232988E+00 0.308E+04 -0.37557232988E+00 0.00E+00  epot                
  0.76705444896E-01 0.308E+04  0.76705444896E-01 0.00E+00  ekin                
  0.89398908362E-01 0.308E+04  0.89398908362E-01 0.00E+00  e2                  
 ===>> dmc block            2
 -0.30052958711E+00 0.308E+04 -0.29969863985E+00 0.83E-03  elocal              
  0.97127499268E+00 0.308E+04  0.97203960103E+00 0.76E-03  acc.rate            
 -0.37623738972E+00 0.308E+04 -0.37590502132E+00 0.33E-03  epot                
  0.75707802610E-01 0.308E+04  0.76206381462E-01 0.50E-03  ekin                
  0.90482170036E-01 0.308E+04  0.89940802283E-01 0.54E-03  e2                  
 ===>> dmc block            3
 -0.30000136817E+00 0.303E+04 -0.29979839419E+00 0.49E-03  elocal              
  0.97175968465E+00 0.303E+04  0.97194736364E+00 0.45E-03  acc.rate            
 -0.37486189202E+00 0.303E+04 -0.37556129177E+00 0.40E-03  epot                
  0.74860523852E-01 0.303E+04  0.75762897585E-01 0.53E-03  ekin                
  0.90111658307E-01 0.303E+04  0.89997102364E-01 0.32E-03  e2                  
 ===>> dmc block            4
 -0.29918209523E+00 0.299E+04 -0.29964695477E+00 0.38E-03  elocal              
  0.96800531437E+00 0.299E+04  0.97097870773E+00 0.10E-02  acc.rate            
 -0.37636514392E+00 0.299E+04 -0.37575881750E+00 0.34E-03  epot                
  0.77183048695E-01 0.299E+04  0.76111862724E-01 0.52E-03  ekin                
  0.89597866902E-01 0.299E+04  0.89899000649E-01 0.25E-03  e2                  
 ===>> dmc block            5
 -0.29963629924E+00 0.302E+04 -0.29964483853E+00 0.30E-03  elocal              
  0.97374771447E+00 0.302E+04  0.97152864612E+00 0.97E-03  acc.rate            
 -0.37732445215E+00 0.302E+04 -0.37606976026E+00 0.41E-03  epot                
  0.77688152910E-01 0.302E+04  0.76424921727E-01 0.51E-03  ekin                
  0.89856292880E-01 0.302E+04  0.89890518676E-01 0.19E-03  e2                  
 ===>> dmc block            6
 -0.29898611972E+00 0.299E+04 -0.29953654318E+00 0.26E-03  elocal              
  0.97470202508E+00 0.299E+04  0.97205035915E+00 0.95E-03  acc.rate            
 -0.37655878569E+00 0.299E+04 -0.37615015751E+00 0.35E-03  epot                
  0.77572665966E-01 0.299E+04  0.76613614328E-01 0.46E-03  ekin                
  0.89460226313E-01 0.299E+04  0.89819777319E-01 0.17E-03  e2                  
 ===>> dmc block            7
 -0.29839705900E+00 0.296E+04 -0.29937683916E+00 0.28E-03  elocal              
  0.97469357114E+00 0.296E+04  0.97242081763E+00 0.89E-03  acc.rate            
 -0.37458684577E+00 0.296E+04 -0.37593105208E+00 0.37E-03  epot                
  0.76189786771E-01 0.296E+04  0.76554212920E-01 0.39E-03  ekin                
  0.89108528260E-01 0.296E+04  0.89720092448E-01 0.18E-03  e2                  
 ===>> dmc block            8
 -0.29897502745E+00 0.300E+04 -0.29932694783E+00 0.24E-03  elocal              
  0.97379842309E+00 0.300E+04  0.97259186932E+00 0.79E-03  acc.rate            
 -0.37158713110E+00 0.300E+04 -0.37539168501E+00 0.63E-03  epot                
  0.72612103653E-01 0.300E+04  0.76064737185E-01 0.60E-03  ekin                
  0.89467382814E-01 0.300E+04  0.89688714517E-01 0.16E-03  e2                  
 ===>> dmc block            9
 -0.30031588700E+00 0.305E+04 -0.29943788161E+00 0.24E-03  elocal              
  0.96933382028E+00 0.305E+04  0.97222639922E+00 0.78E-03  acc.rate            
 -0.37460230707E+00 0.305E+04 -0.37530313692E+00 0.56E-03  epot                
  0.74286420069E-01 0.305E+04  0.75865255311E-01 0.56E-03  ekin                
  0.90296511997E-01 0.305E+04  0.89756893909E-01 0.15E-03  e2                  
 ===>> dmc block           10
 -0.29949322009E+00 0.302E+04 -0.29944340686E+00 0.22E-03  elocal              
  0.96910551847E+00 0.302E+04  0.97191479609E+00 0.77E-03  acc.rate            
 -0.37546966024E+00 0.302E+04 -0.37531976338E+00 0.50E-03  epot                
  0.75976440150E-01 0.302E+04  0.75876356519E-01 0.50E-03  ekin                
  0.89771325064E-01 0.302E+04  0.89758334782E-01 0.14E-03  e2                  
 ===>> dmc block           11
 -0.29860860712E+00 0.298E+04 -0.29936853017E+00 0.21E-03  elocal              
  0.97399101862E+00 0.298E+04  0.97210102120E+00 0.72E-03  acc.rate            
 -0.37524306565E+00 0.298E+04 -0.37531288404E+00 0.45E-03  epot                
  0.76634458533E-01 0.298E+04  0.75944353868E-01 0.46E-03  ekin                
  0.89239801752E-01 0.298E+04  0.89711825380E-01 0.13E-03  e2                  
 ===>> dmc block           12
 -0.29661087123E+00 0.290E+04 -0.29914713866E+00 0.30E-03  elocal              
  0.97575082431E+00 0.290E+04  0.97239403619E+00 0.72E-03  acc.rate            
 -0.37551374797E+00 0.290E+04 -0.37532900988E+00 0.42E-03  epot                
  0.78902876740E-01 0.290E+04  0.76181871213E-01 0.49E-03  ekin                
  0.88031486440E-01 0.290E+04  0.89576923721E-01 0.18E-03  e2                  
 ===>> dmc block           13
 -0.29800787995E+00 0.296E+04 -0.29906087382E+00 0.29E-03  elocal              
  0.97250300220E+00 0.296E+04  0.97240228711E+00 0.67E-03  acc.rate            
 -0.37485070036E+00 0.296E+04 -0.37529279221E+00 0.38E-03  epot                
  0.76842820413E-01 0.296E+04  0.76231918384E-01 0.45E-03  ekin                
  0.88890789822E-01 0.296E+04  0.89524969560E-01 0.18E-03  e2                  
 ===>> dmc block           14
 -0.29895230449E+00 0.301E+04 -0.29905311360E+00 0.27E-03  elocal              
  0.97332417912E+00 0.301E+04  0.97246818133E+00 0.62E-03  acc.rate            
 -0.37859892811E+00 0.301E+04 -0.37552910540E+00 0.43E-03  epot                
  0.79646623621E-01 0.301E+04  0.76475991800E-01 0.48E-03  ekin                
  0.89453945304E-01 0.301E+04  0.89519892948E-01 0.16E-03  e2                  
 ===>> dmc block           15
 -0.29903131565E+00 0.301E+04 -0.29905165779E+00 0.25E-03  elocal              
  0.97040198994E+00 0.301E+04  0.97233018797E+00 0.59E-03  acc.rate            
 -0.37299018034E+00 0.301E+04 -0.37535953988E+00 0.43E-03  epot                
  0.73958864689E-01 0.301E+04  0.76307882092E-01 0.48E-03  ekin                
  0.89487334083E-01 0.301E+04  0.89517718461E-01 0.15E-03  e2                  
 ===>> dmc block           16
 -0.29931506241E+00 0.303E+04 -0.29906822644E+00 0.23E-03  elocal              
  0.97295032639E+00 0.303E+04  0.97236919584E+00 0.56E-03  acc.rate            
 -0.37682798735E+00 0.303E+04 -0.37545190800E+00 0.41E-03  epot                
  0.77512924940E-01 0.303E+04  0.76383681564E-01 0.46E-03  ekin                
  0.89654185633E-01 0.303E+04  0.89526302504E-01 0.14E-03  e2                  
 ===>> dmc block           17
 -0.29906554483E+00 0.301E+04 -0.29906806851E+00 0.22E-03  elocal              
  0.97210217711E+00 0.301E+04  0.97235347031E+00 0.52E-03  acc.rate            
 -0.37754300833E+00 0.301E+04 -0.37557505918E+00 0.41E-03  epot                
  0.78477463498E-01 0.301E+04  0.76506990668E-01 0.44E-03  ekin                
  0.89527848419E-01 0.301E+04  0.89526393547E-01 0.13E-03  e2                  
 ===>> dmc block           18
 -0.29897233396E+00 0.301E+04 -0.29906274707E+00 0.21E-03  elocal              
  0.96967676478E+00 0.301E+04  0.97220468459E+00 0.51E-03  acc.rate            
 -0.37385881543E+00 0.301E+04 -0.37547966110E+00 0.40E-03  epot                
  0.74886481464E-01 0.301E+04  0.76416914034E-01 0.43E-03  ekin                
  0.89478897641E-01 0.301E+04  0.89523753469E-01 0.13E-03  e2                  
 ===>> dmc block           19
 -0.29893534758E+00 0.301E+04 -0.29905604151E+00 0.19E-03  elocal              
  0.97399888836E+00 0.301E+04  0.97229912094E+00 0.50E-03  acc.rate            
 -0.38058142576E+00 0.301E+04 -0.37574818803E+00 0.46E-03  epot                
  0.81646078174E-01 0.301E+04  0.76692146523E-01 0.49E-03  ekin                
  0.89416109662E-01 0.301E+04  0.89518087731E-01 0.12E-03  e2                  
 ===>> dmc block           20
 -0.29958991879E+00 0.303E+04 -0.29908294325E+00 0.19E-03  elocal              
  0.97218722507E+00 0.303E+04  0.97229348258E+00 0.47E-03  acc.rate            
 -0.37792036879E+00 0.303E+04 -0.37585764288E+00 0.45E-03  epot                
  0.78330450002E-01 0.303E+04  0.76774699621E-01 0.47E-03  ekin                
  0.89848532172E-01 0.303E+04  0.89534738622E-01 0.12E-03  e2                  
 ===>> dmc block           21
 -0.29851257512E+00 0.299E+04 -0.29905595922E+00 0.18E-03  elocal              
  0.97205183896E+00 0.299E+04  0.97228205045E+00 0.45E-03  acc.rate            
 -0.37338825856E+00 0.299E+04 -0.37574081664E+00 0.44E-03  epot                
  0.74875683444E-01 0.299E+04  0.76684857423E-01 0.46E-03  ekin                
  0.89187008824E-01 0.299E+04  0.89518287572E-01 0.11E-03  e2                  
 ===>> dmc block           22
 -0.29949019854E+00 0.303E+04 -0.29907583113E+00 0.17E-03  elocal              
  0.97202816235E+00 0.303E+04  0.97227043187E+00 0.43E-03  acc.rate            
 -0.37790708399E+00 0.303E+04 -0.37583995064E+00 0.44E-03  epot                
  0.78416885444E-01 0.303E+04  0.76764119507E-01 0.44E-03  ekin                
  0.89773251672E-01 0.303E+04  0.89529955389E-01 0.11E-03  e2                  
 ===>> dmc block           23
 -0.29845939549E+00 0.299E+04 -0.29904918453E+00 0.17E-03  elocal              
  0.97323990382E+00 0.299E+04  0.97231233914E+00 0.41E-03  acc.rate            
 -0.37476031458E+00 0.299E+04 -0.37579328132E+00 0.42E-03  epot                
  0.76300919088E-01 0.299E+04  0.76744096788E-01 0.42E-03  ekin                
  0.89140543236E-01 0.299E+04  0.89513122308E-01 0.10E-03  e2                  
 ===>> dmc block           24
 -0.29942449632E+00 0.303E+04 -0.29906491647E+00 0.16E-03  elocal              
  0.97299678843E+00 0.303E+04  0.97234102920E+00 0.39E-03  acc.rate            
 -0.37544524201E+00 0.303E+04 -0.37577869255E+00 0.40E-03  epot                
  0.76020745684E-01 0.303E+04  0.76713776078E-01 0.41E-03  ekin                
  0.89726013825E-01 0.303E+04  0.89522046083E-01 0.99E-04  e2                  
 ===>> dmc block           25
 -0.29756098488E+00 0.295E+04 -0.29900587929E+00 0.16E-03  elocal              
  0.97392573845E+00 0.295E+04  0.97240323733E+00 0.38E-03  acc.rate            
 -0.37793735936E+00 0.295E+04 -0.37586343152E+00 0.39E-03  epot                
  0.80376374473E-01 0.295E+04  0.76857552231E-01 0.42E-03  ekin                
  0.88599804640E-01 0.295E+04  0.89485843280E-01 0.10E-03  e2                  
 ===>> dmc block           26
 -0.29955376477E+00 0.303E+04 -0.29902711347E+00 0.16E-03  elocal              
  0.97124038309E+00 0.303E+04  0.97235816904E+00 0.37E-03  acc.rate            
 -0.37476223095E+00 0.303E+04 -0.37582075271E+00 0.38E-03  epot                
  0.75208466185E-01 0.303E+04  0.76793639241E-01 0.41E-03  ekin                
  0.89822790964E-01 0.303E+04  0.89498902231E-01 0.98E-04  e2                  
 ===>> dmc block           27
 -0.29838723598E+00 0.298E+04 -0.29900360919E+00 0.16E-03  elocal              
  0.97342242107E+00 0.298E+04  0.97239726165E+00 0.36E-03  acc.rate            
 -0.37319883846E+00 0.298E+04 -0.37572444332E+00 0.38E-03  epot                
  0.74811602486E-01 0.298E+04  0.76720834127E-01 0.40E-03  ekin                
  0.89097480779E-01 0.298E+04  0.89484157028E-01 0.96E-04  e2                  
 ===>> dmc block           28
 -0.29819723911E+00 0.298E+04 -0.29897503565E+00 0.15E-03  elocal              
  0.96984346659E+00 0.298E+04  0.97230676851E+00 0.36E-03  acc.rate            
 -0.37213997997E+00 0.298E+04 -0.37559742868E+00 0.39E-03  epot                
  0.73942740856E-01 0.298E+04  0.76622393031E-01 0.40E-03  ekin                
  0.88997215208E-01 0.298E+04  0.89466902357E-01 0.94E-04  e2                  
 ===>> dmc block           29
 -0.29676812355E+00 0.292E+04 -0.29890107261E+00 0.17E-03  elocal              
  0.97203373129E+00 0.292E+04  0.97229761787E+00 0.35E-03  acc.rate            
 -0.37221034922E+00 0.292E+04 -0.37548391320E+00 0.39E-03  epot                
  0.75442225666E-01 0.292E+04  0.76582840588E-01 0.38E-03  ekin                
  0.88143207082E-01 0.292E+04  0.89422539681E-01 0.10E-03  e2                  
 ===>> dmc block           30
 -0.29852269273E+00 0.299E+04 -0.29888848899E+00 0.16E-03  elocal              
  0.97547803440E+00 0.299E+04  0.97240338763E+00 0.35E-03  acc.rate            
 -0.37258970395E+00 0.299E+04 -0.37538766172E+00 0.39E-03  epot                
  0.74067011219E-01 0.299E+04  0.76499172721E-01 0.38E-03  ekin                
  0.89216541395E-01 0.299E+04  0.89415688884E-01 0.98E-04  e2                  
 ===>> dmc block           31
 -0.29847236141E+00 0.299E+04 -0.29887511684E+00 0.16E-03  elocal              
  0.97203956984E+00 0.299E+04  0.97239169643E+00 0.34E-03  acc.rate            
 -0.37308918539E+00 0.299E+04 -0.37531380073E+00 0.38E-03  epot                
  0.74616823981E-01 0.299E+04  0.76438683898E-01 0.37E-03  ekin                
  0.89177612441E-01 0.299E+04  0.89408038355E-01 0.95E-04  e2                  
 ===>> dmc block           32
 -0.29886956184E+00 0.300E+04 -0.29887494328E+00 0.15E-03  elocal              
  0.97045928322E+00 0.300E+04  0.97233132216E+00 0.33E-03  acc.rate            
 -0.37542501074E+00 0.300E+04 -0.37531727526E+00 0.37E-03  epot                
  0.76555448905E-01 0.300E+04  0.76442331980E-01 0.36E-03  ekin                
  0.89421091407E-01 0.300E+04  0.89408446171E-01 0.92E-04  e2                  
 ===>> dmc block           33
 -0.29845870558E+00 0.299E+04 -0.29886235420E+00 0.15E-03  elocal              
  0.96946319730E+00 0.299E+04  0.97224457592E+00 0.33E-03  acc.rate            
 -0.37444702903E+00 0.299E+04 -0.37529095472E+00 0.36E-03  epot                
  0.75988323455E-01 0.299E+04  0.76428600522E-01 0.35E-03  ekin                
  0.89147446409E-01 0.299E+04  0.89400552250E-01 0.90E-04  e2                  
 ===>> dmc block           34
 -0.29990114962E+00 0.305E+04 -0.29889339152E+00 0.15E-03  elocal              
  0.96947524928E+00 0.305E+04  0.97216183346E+00 0.33E-03  acc.rate            
 -0.37453695542E+00 0.305E+04 -0.37526842659E+00 0.35E-03  epot                
  0.74635805804E-01 0.305E+04  0.76375035068E-01 0.34E-03  ekin                
  0.90063419288E-01 0.305E+04  0.89420357515E-01 0.89E-04  e2                  
 ===>> dmc block           35
 -0.29820437291E+00 0.298E+04 -0.29887381823E+00 0.14E-03  elocal              
  0.97397280966E+00 0.298E+04  0.97221327875E+00 0.33E-03  acc.rate            
 -0.37630131256E+00 0.298E+04 -0.37529776829E+00 0.34E-03  epot                
  0.78096939654E-01 0.298E+04  0.76423950058E-01 0.34E-03  ekin                
  0.88994388772E-01 0.298E+04  0.89408256811E-01 0.88E-04  e2                  
 ===>> dmc block           36
 -0.29737318281E+00 0.294E+04 -0.29883299450E+00 0.14E-03  elocal              
  0.96985540402E+00 0.294E+04  0.97214913443E+00 0.33E-03  acc.rate            
 -0.37494580202E+00 0.294E+04 -0.37528819330E+00 0.33E-03  epot                
  0.77572619216E-01 0.294E+04  0.76455198793E-01 0.33E-03  ekin                
  0.88499732558E-01 0.294E+04  0.89383541049E-01 0.89E-04  e2                  
 ===>> dmc block           37
 -0.29802164074E+00 0.297E+04 -0.29881125859E+00 0.14E-03  elocal              
  0.97056360104E+00 0.297E+04  0.97210665849E+00 0.32E-03  acc.rate            
 -0.37353659000E+00 0.297E+04 -0.37524126839E+00 0.33E-03  epot                
  0.75514949254E-01 0.297E+04  0.76430009802E-01 0.32E-03  ekin                
  0.88912836930E-01 0.297E+04  0.89370931033E-01 0.87E-04  e2                  
 ===>> dmc block           38
 -0.29914619138E+00 0.303E+04 -0.29882014825E+00 0.14E-03  elocal              
  0.97182850293E+00 0.303E+04  0.97209927578E+00 0.31E-03  acc.rate            
 -0.37512647072E+00 0.303E+04 -0.37523822147E+00 0.32E-03  epot                
  0.75980279338E-01 0.303E+04  0.76418073219E-01 0.31E-03  ekin                
  0.89618180581E-01 0.303E+04  0.89377493441E-01 0.85E-04  e2                  
 ===>> dmc block           39
 -0.29926667225E+00 0.303E+04 -0.29883169702E+00 0.14E-03  elocal              
  0.97364580829E+00 0.303E+04  0.97213927486E+00 0.31E-03  acc.rate            
 -0.37491075767E+00 0.303E+04 -0.37522975204E+00 0.31E-03  epot                
  0.75644085421E-01 0.303E+04  0.76398055018E-01 0.31E-03  ekin                
  0.89632395315E-01 0.303E+04  0.89384086151E-01 0.83E-04  e2                  
 ===>> dmc block           40
 -0.29985404690E+00 0.307E+04 -0.29885780140E+00 0.13E-03  elocal              
  0.97375659507E+00 0.307E+04  0.97218057103E+00 0.30E-03  acc.rate            
 -0.37768242290E+00 0.307E+04 -0.37529237780E+00 0.31E-03  epot                
  0.77828375996E-01 0.307E+04  0.76434576401E-01 0.30E-03  ekin                
  0.89972382608E-01 0.307E+04  0.89399107534E-01 0.82E-04  e2                  
  max.multiplicity           20  exceeded            0  times