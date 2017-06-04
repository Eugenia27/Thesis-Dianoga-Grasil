function looking_progenitors, ireg, iflav, isnapini, isnapfin, onestep_flag, dobh_flag, idhalo, wholecluster_flag, stars_only_flag, misub, rmax_starsonly, yhist, mass
  ; Returns an array with the main progenitor of a given OBJECT at each snap between isnapini and isnapfin, isnapini > isnapfin, if onestep_flag is set then only at the mentioned snaps.
  ; Follows back in time a set of ids particles serching for progenitor of cluster idhalo, or subhalo "idhalo+misub" or "id sub" as given by subfind,
  ; if idhalo < 0 then misub must to be understood as the indexx in the blocks of subfind (idhalo can be any negative number, it is ignored since the whole information is in misub)
  ; if idhalo >= 0 then misub is relative to the first sub of each halo, fsub. In this case you must specify idhalo
  ; onestep_flag        if not 0, then only one step is made between isnapini and isnapfin
  ; wholecluster_flag   if not 0, then all particles in the halo will be tracked (regardless misub), ergo idhalo must be positive, if 0 then particles in misub are tracked, used in ids()
  ; ireg                number of region, from 0 to 28
  
  ;stars_only_flag = 0          ; if defined (not 0) then only stars particles are followed back
  ;rmax_starsonly  = -1        ; radius inside which stars are selected to follow them back, if -1 then all stars in the object are to be tracked
  
  mass = fltarr(isnapini-isnapfin+1)  ; stellar mass inside rmax_starsonly

  path_local     = '/pico/scratch/userexternal/cragonef/tests/'
  path_snapshots = ['/gss/gss_work/DRES_murante/CLUSTERS/PassiveTracers/']
  path_snapshots = ['/gss/gss_work/DRES_murante/CLUSTERS/Dianoga/']
  path_bhdata    = ['/pico/scratch/userexternal/cragonef/tests/']
  
  regiones=['none','D1','D2','D3','D4','D5','D6','D7','D8','D9','D10',     $
           'D11','D12','D13','D14','D15','D16','D17','D18','D19','D20',    $
           'D21','D22','D23','D24','D25','D26','D27','D28','D29']
  flavs=['BH2015']
  
  print, '--------------- SEARCHING PROGENITOR ----------------'
  nomrad = '/pico/scratch/userexternal/cragonef/mferraro/g3d/DirectProgenitor/'
  nomfil = nomrad+'DirectProgenitor_regD'+strtrim(string(ireg),2)+'.dat' 
  openw, unout, nomfil, /GET_LUN
  FM     = '(2i8,3x,3f11.3,3E12.3)' 
  FMdobh = '(2i8,3x,3f11.3,3E12.3,g12.4, f9.3, 1x, i12)'
  
  jj = iflav
  
  ;--- alguna inicializacion se hace fuera del lazo for, se ven los files a isnapini, ie OBJECT que se quiere seguir
  subfile     = path_snapshots+regiones[ireg]+'/'+flavs[jj]+'/Subfind/groups_'+snap_int_to_char(isnapini)+'/sub_'+snap_int_to_char(isnapini)
  readnew, subfile, nsub   , 'NSUB'   ; num de subhalos por halo
  readnew, subfile, fsub   , 'FSUB'   ; first subhalo de cada halo
  readnew, subfile, grnr   , 'GRNR'   ; grupo al que pertenece un dado subhalo
  readnew, subfile, soff   , 'SOFF'   ; indice (no confundir con id) en pid de las particulas en subhalos, es siempre creciente
  readnew, subfile, spos   , 'SPOS'   ; posiciones de los subhalos
  readnew, subfile, smst   , 'SMST'   ; [6,nsubtot] masa de cada tipo para cada subhalo
  readnew, subfile, ncon   , 'NCON'   ; Si es distinto de 0 el HALO esta contaminado por particulas de la zona de baja resolucion
  readnew, subfile, pid    , 'PID '   ; ids de las partículas...
  
  ;--- Inizialize the main_progenitor array with the id of the halo or the subhalo we wish to follow back
  if onestep_flag gt 0 then begin
    main_progenitor = intarr(2)
    delta_snap      = isnapini-isnapfin
  endif else begin
    main_progenitor = intarr(isnapini-isnapfin+1)
    delta_snap      = 1
  endelse
  
  if wholecluster_flag then begin
    if idhalo lt 0 then begin
      stop
    endif
    main_progenitor[0] = idhalo
    print,  '.......', isnapini, main_progenitor[0],  spos[*,fsub[main_progenitor[0]]]
    if(ncon[idhalo] ne 0) then print, 'HALO CONTAMINADO'
  endif else begin
    if idhalo ge 0 then begin
      main_progenitor[0] = fsub[idhalo]+misub
      if(ncon[idhalo] ne 0) then print, 'HALO CONTAMINADO'
    endif else begin
      main_progenitor[0] = misub
      if(ncon[grnr[misub]] ne 0) then print, 'HALO CONTAMINADO'
    endelse
    print, '.......', isnapini, main_progenitor[0],  spos[*,main_progenitor[0]] 
  endelse
  progenitor_clus = idhalo
  
  ;--- Aqui comienza el lazo que  recorre los snaps hacia atras
  ;--- Particulas que pertenecen al main progenitor, seran rastreadas hacia atras
  submass_in_BH = 0.
  for isnap= isnapini, isnapfin+1, -delta_snap do begin
    
    ii = isnapini-isnap
      
    ;--- se selecionan las particulas en el progenitor (o el objeto a rastrear en el primer paso) que da subfind
    if wholecluster_flag then begin
      ids_main_progenitor = particles_ids(main_progenitor[ii], wholecluster_flag, -999, soff, fsub, pid)
      pos = spos[*,fsub[main_progenitor[ii]]]
      ;--- si el subhalo principal tiene asociada masaBH sobreescribe submass_in_BH
      if dobh_flag then submass_in_BH = smst[5,fsub[main_progenitor[ii]]]
    endif else begin
      ; if idhalo < 0 (-999 in this case) then the third field of ids (mi sub) must to be understood as the indexx in the blocks of subfind,
      ; otherwise  is relative to the first sub of each halo, fsub
      ids_main_progenitor = particles_ids(-999, wholecluster_flag, main_progenitor[ii], soff, fsub, pid)
      pos = spos[*,main_progenitor[ii]]
      ;--- si el subhalo tiene asociada masaBH sobreescribe submass_in_BH
      if dobh_flag then submass_in_BH = smst[5,main_progenitor[ii]]
    endelse
    
    if stars_only_flag then begin
      ;--- Se quiere rastrear solo estrellas. Cuales particulas de ids_main_progenitor son estrellas?
      ;--- Usando el snapshot, primero separamos del ids_main_progenitor dado, las particulas que son estrellas
      nombre_snap = path_snapshots+regiones[ireg]+'/'+flavs[jj]+'/snap_'+snap_int_to_char(isnap)
      readnew, nombre_snap, head     ,'HEAD'
      readnew, nombre_snap, ids_stars,'ID  ', parttype= 4    
      readnew, nombre_snap, posstars ,'POS ', parttype= 4
      readnew, nombre_snap, mass_s   ,'MASS', parttype= 4, /QUIET
      if dobh_flag then begin
        readnew, nombre_snap, bhma     ,'BHMA'   ; theoretical bh mass from accretion rates, not dynamical, which is MASS (en flavs NoSwal son la misma)
        readnew, nombre_snap, posbhs   ,'POS ', parttype= 5
        readnew, nombre_snap, ids_bhs  ,'ID  ', parttype= 5
      endif
      
      if rmax_starsonly gt 0 then begin
        dstars = sqrt((posstars[0,*]-pos[0])^2 + (posstars[1,*]-pos[1])^2 + (posstars[2,*]-pos[2])^2)
        inds_stars_dentro = where(dstars le rmax_starsonly, Ndentro)
        mass[ii] = total(mass_s[inds_stars_dentro])*1e10/0.72
        ; find_coincident_ids devuelve en inds_myobj los indices del vector ids_stars en donde ids_stars=ids_main_progenitor
        ;--- si solo las estrellas dentro de rmax_starsonly
        FIND_COINCIDENT_IDS, ids_stars[inds_stars_dentro], ids_main_progenitor, inds_myobj   
      endif else begin
        ;--- si todas las estrellas en el subhalo
        FIND_COINCIDENT_IDS, ids_stars, ids_main_progenitor, inds_myobj
      endelse
    endif
    
    if ii eq 0 then begin 
      if dobh_flag  then begin
        idBH_in_sub   = 0
        dfirstBH      = -1. 
        if submass_in_BH gt 0 then begin
          ;--- Si el subhalo tiene asociado BHs entonces vemos el BH mas cercano
          dbhs     = sqrt((posbhs[0,*]-pos[0])^2 + (posbhs[1,*]-pos[1])^2 + (posbhs[2,*]-pos[2])^2)
          indord   = sort(dbhs)
          dfirstBH = dbhs[indord[0]]
          ;--- reescribo submass_in_BH con la masa del BH mas cercano al subhalo
          submass_in_BH = bhma[indord[0]]*1e10
          idBH_in_sub   = ids_bhs[indord[0]]
        endif
      endif else begin
      endelse    
      flush , unout
    endif

    ;--- Lectura z anterior, Buscamos a este nuevo z que posicion en pid tienen las particulas ids_main_progenitor del z anterior
    isnapprev= isnap-delta_snap
    subfile = path_snapshots+regiones[ireg]+'/'+flavs[jj]+'/Subfind/groups_'+snap_int_to_char(isnapprev)+'/sub_'+snap_int_to_char(isnapprev)
    readnew, subfile, nsub   , 'NSUB'   ; num de subhalos por halo
    readnew, subfile, fsub   , 'FSUB'   ; first subhalo de cada halo
    readnew, subfile, soff   , 'SOFF'   ; indice inicial (no confundir con id) en pid de las particulas en subhalos, es siempre creciente
    readnew, subfile, pid    , 'PID '   ; ids de las partículas...
    readnew, subfile, spos   , 'SPOS'
    readnew, subfile, grnr   , 'GRNR'
    readnew, subfile, smst   , 'SMST'
    nombre_snap = path_snapshots+regiones[ireg]+'/'+flavs[jj]+'/snap_'+snap_int_to_char(isnapprev)
    readnew, nombre_snap, head     ,'HEAD'
    if dobh_flag then begin
      readnew, nombre_snap, bhma     ,'BHMA'   ; theoretical bh mass from accretion rates, not dynamical, which is MASS (en flavs NoSwal son la misma)
      readnew, nombre_snap, posbhs   ,'POS ', parttype= 5
      readnew, nombre_snap, ids_bhs  ,'ID  ', parttype= 5
    endif
    ;--- busco los indices de pid que contienen a las particulas del paso anterior
    if stars_only_flag then begin
      if rmax_starsonly gt 0 then begin
        ;--- si solo las estrellas dentro de rmax_starsonly
        FIND_COINCIDENT_IDS, pid, ids_stars[inds_stars_dentro[inds_myobj]], inds_en_pid
      endif else begin
        ;--- si todas las estrellas
        FIND_COINCIDENT_IDS, pid, ids_stars[inds_myobj], inds_en_pid
      endelse
    endif else begin
      FIND_COINCIDENT_IDS, pid, ids_main_progenitor, inds_en_pid
    endelse
    
    part_total_x_halo = ulonarr(n_elements(nsub))
    
    ;--- busco finalmente en que subhalo esta cada particula
    inds_en_soff = value_locate(soff,inds_en_pid)
    
    ;--- cuento cuantas particulas caen en cada subhalo, xhist es el subhalo, no todo xhist tiene particulas, i.e, existen yhist=0
    yhist = histogram(inds_en_soff, LOCATIONS= xhist)
    ;plot, xhist, yhist, yrange=[0,100],color=100
    
    ;--- sumo las particulas de los subhalos que pertenecen a un halo
    for ihalo=0 , n_elements(nsub)-1 do begin
      ind = where(grnr[xhist] eq ihalo, nxhalo)
      if nxhalo gt 0 then part_total_x_halo[ihalo] = total(yhist[ind])
    endfor
    
    
    if wholecluster_flag then begin
      ;--- el main_progenitor es el cumulo que contiene mas particulas
      part_total_x_halo_max = max(part_total_x_halo,indmax)
      main_progenitor[ii+1] = indmax
      pos = spos[*,fsub[main_progenitor[ii+1]]]
      print, '.......', isnapprev, main_progenitor[ii+1],  pos
    endif else begin
      ;--- main progenitor es el subhalo que tiene mas particulas
      yhistmax = max(yhist,indmax)
      main_progenitor[ii+1] =  xhist[indmax]
      pos = spos[*,xhist[indmax]]
      ;--- en que cumulo esta este main_progenitor
      progenitor_clus = grnr[xhist[indmax]]
      print, '.......', isnapprev, progenitor_clus, main_progenitor[ii+1],  pos
    endelse
    
    ind_progenitors = where(yhist ne 0, n_progenitors)
    
    for i=0 , n_progenitors-1 do begin  
      if dobh_flag then begin   
        submass_in_BH = smst[5,xhist[ind_progenitors[i]]]
        idBH_in_sub   = 0
        dfirstBH      = -1.        
        if submass_in_BH gt 0 then begin
          ;--- Si el subhalo tiene asociado BHs entonces vemos el BH mas cercano
          dbhs     = sqrt((posbhs[0,*]-spos[0,xhist[ind_progenitors[i]]])^2 + (posbhs[1,*]-spos[1,xhist[ind_progenitors[i]]])^2 + (posbhs[2,*]-spos[2,xhist[ind_progenitors[i]]])^2)
          indord   = sort(dbhs)
          dfirstBH = dbhs[indord[0]]
          ;--- reescribo submass_in_BH con la masa del BH mas cercano al subhalo
          submass_in_BH = bhma[indord[0]]*1e10
          idBH_in_sub   = ids_bhs[indord[0]]
        endif 
      endif else begin
      endelse
      flush , unout
    endfor
    if onestep_flag then break ; will force to leave the for loop
    
  endfor  ; de los snaps
  print, ''
  print, '---------------- END SEARCHING PROGENITOR -----------------'
  print, 'Main Progenitor= ', main_progenitor
  
  
  return, main_progenitor
end

