function particles_ids, idhalo, todoelhalo, misub, soff, fsub, pid
  ; if todoelhalo ne 0,  escribe los ids de todo el halo y no da bola a misub
  ; if idhalo < 0 then misub must to be understood as the indexx in the blocks of subfind, otherwise misub is relative to the first sub of each halo, fsub
  path_local     = '/pico/scratch/userexternal/cragonef/tests/'
  path_snapshots = ['/gss/gss_work/DRES_murante/CLUSTERS/PassiveTracers/']
  path_snapshots = ['/gss/gss_work/DRES_murante/CLUSTERS/Dianoga/']
  path_bhdata    = ['/pico/scratch/userexternal/cragonef/tests/']

  c = 0UL
  print_flag = 0
  
  if print_flag then begin
    file_para_gadgetviewer = '/bigcin/dianoga/ids.dat'
    openw, 33, file_para_gadgetviewer
    print, file_para_gadgetviewer
  endif
  
  
  if todoelhalo then begin
    ;indp_sub  = fsub[idhalo]                                    ; indice de la primer sub perteneciente a indhalo
    ;indu_sub  = fsub[idhalo+1]-1                                ; indice de la ultima sub perteneciente a indhalo
    ;ini       = soff[indp_sub]                                  ; SI SE QUIEREN TODAS LAS PARTICULAS DEL HALO
    ;fin       = soff[indu_sub]-1                                ; SI SE QUIEREN TODAS LAS PARTICULAS DEL HALO
    
    ini       = soff[fsub[idhalo]]
    fin       = soff[fsub[idhalo+1]]-1
  endif else begin
    if idhalo lt 0 then begin
      misubl = misub
    endif else begin
      misubl = fsub[idhalo]+misub
    endelse 
    ini       = soff[misubl]
    fin       = soff[misubl+1]-1
  endelse
  
  my_ids = ulonarr(fin-ini+1)
  
  for i= ini, fin do begin
    if print_flag then printf,33, pid[i]
    my_ids[c] = pid[i]
    c++
  endfor
  print, 'Num de IDs de todas las particulas', c
  if print_flag then close, 33
  return, my_ids
end
