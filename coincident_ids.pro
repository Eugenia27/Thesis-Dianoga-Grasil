pro coincident_ids, ids_stars, ids_misub, inds_misub_en_stars, status
  ; used by followback_ids
  
  inds_misub_en_stars = ulonarr(n_elements(ids_misub))
  
  indsord_stars = sort(ids_stars)
  ind_tmp = value_locate(ids_stars[indsord_stars],ids_misub)
  ; dado un ids_misub, value_locate indica entre cuales ids_stars (ordenado) esta, ie: da un intervalo.
  ; Por lo que podria suceder que ids_misub no sea realmente uno de los ids_stars. Por ello se hace un control.
  estoy = 0UL
  for icontrol=0UL, n_elements(ids_misub)-1 do begin
    if ids_misub[icontrol] eq ids_stars[indsord_stars[ind_tmp[icontrol]]] then begin
      inds_misub_en_stars[estoy] = indsord_stars[ind_tmp[icontrol]]
      estoy++
    endif
  endfor
  
  if estoy gt 0 then begin
    ; estos son los indices del vector ids_stars en donde se encuentran las estrellas que pertenecen a misub
    inds_misub_en_stars = inds_misub_en_stars[0:estoy-1]
    status = 1  
  endif else begin
    print, 'Inside FIND_COINCIDENT_IDS: NO MATCHES FOUND'
    status = 0
  endelse

end
