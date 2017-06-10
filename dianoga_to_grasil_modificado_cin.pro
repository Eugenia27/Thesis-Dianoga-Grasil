pro dia_to_g3d, win, ireg, biased_progenitor, direct_progenitor, masks, write=write

  n_s       = 0L
  n_g       = 0L
  n_bh      = 0L
  box       = 500.     
  go_on_BHs = 1  ; en el caso que se use la mascara y no se encuentren bhs go_on_BHs sera 0
  do_BH     = 1  
  date      = SYSTIME()

  simu_dir   = '/gss/gss_work/DRES_murante/CLUSTERS/Dianoga/'
  output_dir = '/pico/scratch/userexternal/cragonef/mferraro/g3d/Input/' 
  plots_dir  = '/pico/scratch/userexternal/cragonef/mferraro/g3d/plots/' 
  ;simu_dir   = '/bigcin/dianoga/'
  ;output_dir = '/pico'
  print, 'Simu DIR: ', simu_dir

  if ireg EQ 0 then begin
    print, 'D0 DOES NOT EXIST'
    exit
  endif
 
  regions=['none','D1','D2','D3','D4','D5','D6','D7','D8','D9','D10',      $
           'D11','D12','D13','D14','D15','D16','D17','D18','D19','D20',    $
           'D21','D22','D23','D24','D25','D26','D27','D28','D29']

  flav='CSF2015'
  if do_BH then flav='BH2015'
  snap='041'
  
  snap_name = simu_dir+regions[ireg]+'/'+flav+'/snap_'+snap
  sub_name  = simu_dir+regions[ireg]+'/'+flav+'/Subfind/groups_'+snap+'/sub_'+snap
  print, 'SNAP Name: ', snap_name
  print, 'SUB  Name: ', sub_name

  ; Cosmological SNAP data
  readnew, snap_name, head, "HEAD"
  z           = head.redshift
  hubble      = head.hubbleparam
  omega_m     = head.omega0
  omega_lambda= head.omegalambda
  print, 'REDSHIFT: ', z
 
  ; gas
  readnew, snap_name, cldx_g , "CLDX" , parttype=0
  readnew, snap_name, pos_g  , "POS " , parttype=0
  readnew, snap_name, M_g    , "MASS" , parttype=0
  readnew, snap_name, Z_g    , "Zs  " , parttype=0
  readnew, snap_name, rho_g  , "RHO " , parttype=0  
  readnew, snap_name, T_g    , "TEMP" , parttype=0
  readnew, snap_name, idga   , "ID  " , parttype=0
  
  ; stars 
  readnew, snap_name, pos_s  , "POS " , parttype=4
  readnew, snap_name, M_s    , "MASS" , parttype=4
  readnew, snap_name, Mi_s   , "iM  " , parttype=4
  readnew, snap_name, Z_s    , "Zs  " , parttype=4
  readnew, snap_name, age_s  , "AGE " , parttype=4
  readnew, snap_name, idst   , "ID  " , parttype=4
  
  ; BHs
  if do_BH then begin
    print, 'AGN-ON'
    readnew, snap_name, pos_bh , "POS " , parttype=5
    readnew, snap_name, M_bh   , "BHMA" , parttype=5
    readnew, snap_name, Mdot_bh, "BHMD ", parttype=5
    readnew, snap_name, idbh   , "ID  " , parttype=5
  endif
  
 
  ; SubFind bolcks 
  readnew, sub_name, spos, 'SPOS'  ; posizioni di tutti i sottoaloni
  readnew, sub_name, gpos, 'GPOS'  ; posizioni di tutti i sottoaloni
  readnew, sub_name, smst, 'SMST'  ; masas de cada subhalo
  readnew, sub_name, fsub, 'FSUB'  ; indici del primo sottoalone di ogni alone negli array di sottoaloni (tipo spos)
  readnew, sub_name, ncon, 'NCON'  ; contaminacion de cada halo, si es distinto de cero esta contaminado 
  readnew, sub_name, grnr, 'GRNR'  ; grupo al que pertenece un dado subhalo
  readnew, sub_name, mtot, 'MTOT'  ; masa total de todas las particulas pertenecientes al subhalo
  readnew, sub_name, m500, 'M500'  ; 
  readnew, sub_name, mvir, 'MVIR'  ;
  readnew, sub_name, msub, 'MSUB'  ; 
  readnew, sub_name, soff, 'SOFF'  ; 
  readnew, sub_name, pid , 'PID '  ; 

  index_no_cont = where(ncon[0:10] eq 0)
  index_max     = where(m500 eq max(m500[index_no_cont]))
  index_group   = index_no_cont[index_max[0]] 
  print, 'MTOT: ', mtot[0:10]   
  print, 'M500: ', m500[0:10]   
  print, 'MVIR: ', mvir[0:10]   
  print, 'MSUB: ', msub[0:10]   
  print, 'NOT CONT GROUP INDEX: ',index_no_cont
  print, 'NOT CONT M500: '       ,m500[index_no_cont]
  print, 'MAX_INDEX:   '         ,index_max
  print, 'GROUP_INDEX: '         ,index_group

  pos_cdm = dblarr(3,2)  
  if biased_progenitor eq 1 then begin
    print,''
    print,'BIASED-PROGENITOR ON'
    b_progenitor = fsub[index_group]
    pos_cdm[*,0] = spos(*,b_progenitor)

    fnameb='biased_progenitors.dat'
    openw,  25, fnameb, /append
    printf, 25, 'D'+strtrim(string(ireg),2), strtrim(string(z,format='(f0.3)'),2),flav, b_progenitor, index_group, fsub[index_group], $
                 date, FORMAT='(A4,2X,A7,2X,A7,2X,I5,2X,I5,2X,I5,2X,A35)'
    close,  25
  endif

  if direct_progenitor eq 1 then begin
    print,''
    print,'DIRECT-PROGENITOR ON'
    progenitors   = looking_progenitors(ireg, 0, 91, UINT(snap), 1, 1, 0, 0, 1, 0, 15, yhist, mass)
    d_progenitor  = progenitors[1]
    group         = grnr[d_progenitor]
    if ncon[group] ne 0 then print, "THE PROGENITOR GROUP IS CONTAMINATED!!!!!!!!!!!!!!!!!!!!!!!!"
    main_of_group = fsub[group]

    fnamed='direct_progenitors.dat'
    openw,  22, fnamed, /append
    printf, 22, 'D'+strtrim(string(ireg),2), strtrim(string(z,format='(f0.3)'),2),flav, d_progenitor, group, main_of_group, $
                 date, FORMAT='(A4,2X,A7,2X,A7,2X,I5,2X,I5,2X,I5,2X,A35)'
    close,  22

    if group eq index_group and d_progenitor eq main_of_group then begin
      print,'DIRECT PROGENITOR = BIASED PROGENITOR; DIRECT PROGENITOR OFF'
      direct_progenitor = 0
    endif else begin
      pos_cdm[*,1] = spos(*,d_progenitor)
      print,'DIRECT PROGENITOR != BIASED PROGENITOR'
    endelse 
  endif
    
  for prog=0,1 do begin
    if prog eq 0 and biased_progenitor eq 0 then continue 
    if prog eq 1 and direct_progenitor eq 0 then continue 

    print,'POS CENTER OF MASS OF PROGENITOR', prog,': ', pos_cdm[*,prog]
    progen_type = ''
    if prog eq 0 then begin
      progen_type = 'BP'
      misub       = b_progenitor
      group_id    = index_group
    endif else begin
      progen_type = 'DP'
      misub       = d_progenitor
      group_id    = group 
    endelse

    box_half_gadget = box*(1+z)*hubble/2.
    file_out_name   = output_dir+flav+'_'+regions[ireg]+'_'+snap+'_'+progen_type+'.dat'    

    ; select only particles in the desired box
    ind_box_g  = where(abs(pos_g[0,*]-pos_cdm[0,prog])  lt box_half_gadget and abs(pos_g[1,*]-pos_cdm[1,prog])  lt box_half_gadget and abs(pos_g[2,*]-pos_cdm[2,prog])  lt box_half_gadget)
    ind_box_s  = where(abs(pos_s[0,*]-pos_cdm[0,prog])  lt box_half_gadget and abs(pos_s[1,*]-pos_cdm[1,prog])  lt box_half_gadget and abs(pos_s[2,*]-pos_cdm[2,prog])  lt box_half_gadget)
    if do_BH then ind_box_bh = where(abs(pos_bh[0,*]-pos_cdm[0,prog]) lt box_half_gadget and abs(pos_bh[1,*]-pos_cdm[1,prog]) lt box_half_gadget and abs(pos_bh[2,*]-pos_cdm[2,prog]) lt box_half_gadget)
    
    if masks eq 0 then begin
      pos_g_  = pos_g[*,ind_box_g]   & M_g_=M_g[ind_box_g]    & Z_g_=Z_g[*,ind_box_g] & rho_g_=rho_g[ind_box_g] & cldx_g_=cldx_g[ind_box_g] & t_g_=t_g[ind_box_g]
      pos_s_  = pos_s[*,ind_box_s]   & M_s_=M_s[ind_box_s]    & Mi_s_=Mi_s[ind_box_s] & Z_s_=Z_s[*,ind_box_s]   & age_s_=age_s[ind_box_s] 
      if do_BH then begin
        pos_bh_ = pos_bh[*,ind_box_bh]
        M_bh_   = M_bh[ind_box_bh] 
        Mdot_bh_= Mdot_bh[ind_box_bh]
      endif 
      name = plots_dir+'D'+strtrim(string(ireg),2)+'_'+flav+'_'+snap+'_'+progen_type+'_XYproj.eps'
    endif else begin
      print, 'INSIDE MASKS '
      idhalo = -1*UINT(grnr(misub))
      ;--- separo todas las particulas que pertenecen a misub
      ids_mask = particles_ids(idhalo, 0, misub, soff, fsub, pid)
      ;--- separo gas, estrellas y bhs. devuelve STATUS=0 si no encuetra matchs 
      COINCIDENT_IDS, idga[ind_box_g] , ids_mask, inds_g_myobj , statusg
      COINCIDENT_IDS, idst[ind_box_s] , ids_mask, inds_s_myobj , status
      if do_BH then begin
        COINCIDENT_IDS, idbh[ind_box_bh], ids_mask, inds_bh_myobj, statusbh
        if statusbh eq 0 then go_on_BHs = 0  ; ya no se consideraran los BHs
      endif
      ;--- finalmente, los indices que estan dentro del box y en misub
      inds_g  = ind_box_g[inds_g_myobj]
      inds_s  = ind_box_s[inds_s_myobj]
      if do_BH and go_on_BHs then inds_bh = ind_box_bh[inds_bh_myobj]
      ;--- redefino vectores
      ; gas
      pos_g_   = pos_g[*,inds_g]
      M_g_     = M_g[inds_g]
      Z_g_     = Z_g[*,inds_g]
      rho_g_   = rho_g[inds_g]
      cldx_g_  = cldx_g[inds_g]
      t_g_     = t_g[inds_g]
      ; estrellas
      pos_s_   = pos_s[*,inds_s]
      M_s_     = M_s[inds_s]
      Mi_s_    = Mi_s[inds_s] 
      Z_s_     = Z_s[*,inds_s]  
      age_s_   = age_s[inds_s] 
      ; bhs
      if do_BH and go_on_BHs then begin
        pos_bh_  = pos_bh[*,inds_bh] 
        M_bh_    = M_bh[inds_bh]
        Mdot_bh_ = Mdot_bh[inds_bh] 
      endif
      name = plots_dir+'D'+strtrim(string(ireg),2)+'_'+flav+'_'+snap+'_'+progen_type+'_XYproj_mask.eps'
    endelse ; fin mascara
  
    ind_dusty = where(cldx_g_ gt 0 or t_g_ lt 1e5, nind_dusty)
    print, 'N DUSTY= ', nind_dusty 
    pos_g_ = pos_g_[*,ind_dusty] & M_g_=M_g_[ind_dusty] & Z_g_=Z_g_[*,ind_dusty] & rho_g_=rho_g_[ind_dusty] & cldx_g_=cldx_g[ind_dusty]
  
    for i=0,2 do begin                   ; refer coordinates to centre of selected object
      pos_g_[i,*]  = pos_g_[i,*] -pos_cdm[i,prog]
      pos_s_[i,*]  = pos_s_[i,*] -pos_cdm[i,prog]
       if do_BH and go_on_BHs then pos_bh_[i,*] = pos_bh_[i,*]-pos_cdm[i,prog]
    endfor
   
    n_g   = n_elements(M_g_)
    n_s   = n_elements(M_s_)
    if do_BH and go_on_BHs then n_bh  = n_elements(M_bh_)
    n_tot = n_g+n_s+n_bh  ; n_bh fue inicializada como cero por lo que aqui no suma si go_on_BHs=0
  
    ; units conversions here only Mdot_bh in order to have the possibility to select BH particles to plot
    if do_BH and go_on_BHs then Mdot_bh_ = Mdot_bh_*1e10/hubble*0.98/1e9  

    if 1 then begin
      dev = 'ps'
      win = 0
      dev = STRUPCASE(dev)
      set_plot, dev
      DEVICE, FILENAME= name, DECOMPOSED=1, /color, /encapsulated,xsize=20, ysize=18, /portrait
      plot,[0],[0], XRANGE= [-250,250], YRANGE= [-250,250], XTITLE= 'x [kpc]', YTITLE= 'y [kpc]',                      $
          TITLE= 'D'+strtrim(string(ireg),2)+' '+progen_type+' z='+strtrim(string(z,format='(f0.3)'),2)+               $
          ' PR='+strtrim(string(misub),2)+' GR='+strtrim(string(group_id),2)+' FS='+strtrim(string(fsub[group_id]),2), $         
          CHARSIZE= 1.7, CHARTHICK= 2, XTHICK= 2, YTHICK= 2, xstyle=1, ystyle=1, /nodata
      oplot,pos_s_[0,*]/(1+z)/hubble, pos_s_[1,*]/(1+z)/hubble, psym=3, color=cgColor('dodger blue')
      oplot,pos_g_[0,*]/(1+z)/hubble, pos_g_[1,*]/(1+z)/hubble, psym=3, color=cgColor('red')
      if do_BH and go_on_BHs then oplot,pos_BH_[0,*]/(1+z)/hubble, pos_BH_[1,*]/(1+z)/hubble, psym=2, color=cgColor('black'), SYMSIZE= 0.7, THICK= 2
      DEVICE, /CLOSE_FILE
    endif
 
    ; all units conversions
    pos_g_  = pos_g_/(1+z)/hubble
    pos_s_  = pos_s_/(1+z)/hubble
    if do_BH and go_on_BHs then pos_bh_ = pos_bh_/(1+z)/hubble
    n_ele = n_elements(z_g_)/n_elements(M_g_)        ; number of chemical elements followed in the simulation
  
    met_g = total(Z_g_[1:n_ele-1,*],1)/M_g_
    met_s = total(Z_s_[1:n_ele-1,*],1)/M_s_
    M_g_   = M_g_*1e10/hubble                          
    M_s_   = M_s_*1e10/hubble
    Mi_s_  = Mi_s_*1e10/hubble
    age_s_ = galage(z,100,H0=hubble*100,lambda0=omega_lambda,omega_m=omega_m)-galage(1/age_s_-1, 100,H0=hubble*100,lambda0=omega_lambda,omega_m=omega_m)
    age_s_ = age_s*1e-9
    min_allowed_age = 1e-4
    jovenes = where(age_s_ lt min_allowed_age,njov)
    if njov gt 0 then age_s_[jovenes] = min_allowed_age
    rho_g_ = rho_g_*(1+z)^3*hubble^2*1e10*1e9        
    if do_BH and go_on_BHs then M_bh_  = M_bh_*1e10/hubble
  
    if keyword_set(write) then begin
      print, 'WRITING FILE ', file_out_name
      get_lun, uni
      openw, uni, file_out_name
      ;printf, uni, -666.666 , n_tot,  1        ;Guardo todas las particulas
      printf, uni, -666.666 , n_tot-n_g,  1     ;No guardo particulas de gas-> engaño para desconsiderar polvo
      for i=0L, n_s-1 do begin
        printf, uni, format='(i10,9g15.5)', -1, pos_s_[*,i], age_s_[i], met_s[i], 0.0, Mi_s_[i], M_s_[i], -666.666
      endfor
      ; No polvo: comentar siguiente blucle for 
      ; for i=0L, n_g-1 do begin     
      ;  printf, uni, format='(i10,9g15.5)', 1, pos_g_[*,i], -666.666, met_g[i], rho_g_[i], M_g_[i], M_g_[i], cldx_g_[i] 
      ; endfor 
      
      if do_BH and go_on_BHs then begin
        for i=0L, n_bh-1 do begin
          printf, uni, format='(i10,9g15.5)', -1, pos_bh_[*,i], 666.666, -666.666, Mdot_bh_[i], M_bh_[i], M_bh_[i], -666.666
        endfor
      endif
      close, uni
      free_lun, uni
    endif
    n_g   = 0L 
    n_s   = 0L
    n_bh  = 0L
    n_tot = 0L
    nind_dusty    = 0
    ind_box_g     = 0
    ind_box_s     = 0
    ind_box_dusty = 0
    pos_g_ = [ ]
    pos_s_ = [ ]
    pos_bh_= [ ]
    M_g_   = [ ]                         
    M_s_   = [ ]
    Z_g_   = [ ]
    T_g_   = [ ]
    Z_s_   = [ ]
    rho_g_ = [ ]
    Mi_s_  = [ ]
    age_s_ = [ ]
    jovenes= [ ]
    njov   = 0 
    n_ele  = 0
    met_g  = 0.
    met_s  = 0.
    if do_BH and go_on_BHs then begin
      ind_box_bh    = 0
      pos_bh_= [ ]
      M_bh_  = [ ]
      Modot_bh_= [ ]
    endif 
  endfor
end


