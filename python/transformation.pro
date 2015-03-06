function acsproc_get_lens_status, status 
  case status of
    0: lens_status = ''
    1: lens_status = 'YES'
    2: lens_status = 'MAYBE'
    3: lens_status = 'NO'  
  endcase  
  return, lens_status     
end


pro transformation

; Run acsproc first to generate the initial metastructure
;spawn, 'idl -e "acsproc"'

run_process = [{title:'REFRESH', statuscode:'Downloaded/Archived'}, $
               {title:'CUTOUT', statuscode:'CUTOUT'}, $
               {title:'LCOORDS', statuscode:'LCOORDS'}, $
               {title:'REGISTER', statuscode:'REGISTER'}, $
               {title:'RECTIFY', statuscode:'RECTIFY'}, $
               {title:'COMBINE', statuscode:'COMBINE'}, $
               {title:'INITIAL MASK', statuscode:'INITIAL MASK'}, $
               {title:'FEATURE/JUNK MASK', statuscode:'FEATURE/JUNK MASK'}, $
               {title:'BSPLINE', statuscode:'BSPLINE'}, $
               {title:'BIZZLE', statuscode:'BIZZLE'}, $
               {title:'PHOTMOD', statuscode:'PHOTMOD'}, $
               {title:'RESET ALL', statuscode:'REDUCED'}]



	target_name=['SHARP0435-1223','SHARP0631+0519']


	image_file=['HE0435_nirc2_n_Kp_6x6.fits','B0631_nirc2_n_H_6x6.fits']


	rms_file=['HE0435_nirc2_n_Kp_6x6_rms.fits','B0631_nirc2_n_H_6x6_rms.fits']

   dataroot=getenv('HST_DATAROOT')+'/18/13000'
   spawn, 'ls '+dataroot+'/*/*flt.fits', allvisits
   pos=strpos(allvisits, '13000')
   visits=strmid(allvisits, pos[0]+6, 2)
   ; update the meta structure
   meta=mrdfits(dataroot+'/.metafile.fits', 1, hdr)
   nvisit=n_elements(meta)
   for i=0L, nvisit-1 do begin &$
      meta[i].TARGNAME=target_name[i] &$
      meta[i].UNIQNAME=target_name[i] &$
      meta[i].VISIT=strtrim(visits[i], 2) &$
      meta[i].PROPOSID='13000' &$
      meta[i].CUTFLAG=1 &$
      meta[i].RECTIFLAG=1 &$
      meta[i].COMBFLAG=1 &$
      maskfile=dataroot+'/'+meta[i].VISIT+'/'+meta[i].TARGNAME+'_13000_55768_F814W_1_mask.fits' &$
      if (file_test(maskfile)) then begin &$
         maskstruc=mrdfits(maskfile, 1, mask_header, /silent) &$
	 photomod_file=dataroot+'/'+meta[i].VISIT+'/'+meta[i].TARGNAME+'_13000_55768_F814W_1_photmod.fits' &$
	 biz_file=dataroot+'/'+meta[i].VISIT+'/'+meta[i].TARGNAME+'_13000_55768_F814W_1_biz.fits' &$
         if maskstruc.has_photmod and file_test(photomod_file) then meta[i].maskflag=4 $
         else if maskstruc.has_bizzle and file_test(biz_file) then meta[i].maskflag=3 $
         else if maskstruc.has_modim then meta[i].maskflag=2 $
         else meta[i].maskflag=1 &$
         void=where(maskstruc.fmask ne 0,fmaskflag) &$
         void=where(maskstruc.jmask ne 0,jmaskflag) &$
         meta[i].fmaskflag=(fmaskflag gt 0) &$
         meta[i].jmaskflag=(jmaskflag gt 0) &$
         meta[i].lens=acsproc_get_lens_status(maskstruc.lens_status) &$
      endif else begin &$
         meta[i].maskflag=file_test(maskfile) &$
         meta[i].fmaskflag=0L &$
         meta[i].jmaskflag=0L &$
         meta[i].lens = acsproc_get_lens_status(0L) &$
      endelse &$

      if (meta[i].maskflag eq 4) then meta[i].statuscode=11 $
      else if (meta[i].maskflag eq 3) then meta[i].statuscode=10 $
      else if (meta[i].maskflag eq 2) then meta[i].statuscode=9 $
      else if (meta[i].fmaskflag) or (meta[i].jmaskflag) then meta[i].statuscode=8  $
      else if (meta[i].maskflag eq 1) then meta[i].statuscode=7 $
      else if (meta[i].combflag) then meta[i].statuscode=6 &$
      meta[i].status=run_process[meta[i].statuscode].statuscode &$
   endfor
   mwrfits, meta, dataroot+'/.metafile.fits', hdr, /create

   for i=0L, nvisit-1 do begin &$
      visit_path=dataroot+'/'+meta[i].VISIT+'/' &$
      ; Update the flt files 
      flt_name=visit_path+'jbjh1chcq_flt.fits' &$
      f=mrdfits(flt_name, 0, hdr_flt, /silent) &$
      sxaddpar, hdr_flt, 'PROPOSID', '13000' &$
      sxaddpar, hdr_flt, 'LINENUM', meta[i].visit+'.001' &$
      sxaddpar, hdr_flt, 'TARGNAME', meta[i].targname &$
      modfits, flt_name, 0, hdr_flt &$

      ; Update the comb file, inserting Keck images
      comb_name=visit_path+target_name[i]+'_13000_55768_F814W_1.fits' &$
      if (NOT file_test(comb_name)) then begin &$
	 print, visits[i] &$
         cmd='cp '+getenv('HST_DIR')+'/18/13000/SLACSJ1541+3642_12210_55734_F814W_1.fits '+comb_name &$


         spawn, cmd &$
         hdr0=headfits(comb_name, ext=0) &$
         old_img=mrdfits(comb_name, 1, hdr1) &$
         old_err=mrdfits(comb_name, 2) &$
         old_mask=mrdfits(comb_name, 3) &$
         old_psf=mrdfits(comb_name, 4, hdr_PSF) &$
         nx=(size(old_img))[1] &$
         ny=(size(old_img))[2] &$

         ; Read in the SHARP target file
         new_img=fltarr(nx, ny) &$
         new_mask=intarr(nx, ny) &$
         new_err=fltarr(nx, ny) &$
         img=mrdfits(visit_path+image_file[i], 0, hdr_target) &$      
         err=mrdfits(visit_path+rms_file[i], 0) &$
         wh=where(img lt 0.0) &$
         img[wh]=0.0 &$
         err[wh]=0.0 &$
         mask=(img gt 0) &$
         xcen=nx/2 & ycen=ny/2 &$
         img_nx=(size(img))[1] & img_ny=(size(img))[2] &$
         new_img[xcen-img_nx/2:xcen+img_nx/2-1, ycen-img_ny/2:ycen+img_ny/2-1]=img &$
         new_err[xcen-img_nx/2:xcen+img_nx/2-1, ycen-img_ny/2:ycen+img_ny/2-1]=err &$
         new_mask[xcen-img_nx/2:xcen+img_nx/2-1, ycen-img_ny/2:ycen+img_ny/2-1]=mask &$
         new_psf=mrdfits(getenv('HST_DIR')+'/18/13000/B1422_PSF.fits', 0) &$
         
         mwrfits, 0, comb_name, hdr0, /create &$
         mwrfits, new_img, comb_name, hdr1 &$
         mwrfits, new_err, comb_name &$
         mwrfits, new_mask, comb_name &$
         mwrfits, new_psf, comb_name &$
      endif &$
   endfor

;   spawn, 'idl -e "acsproc"' ; should be in the "Initial Mask" stage
   		    ; Need to make 13000 the default
end
