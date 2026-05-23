
      itot=itot+1;                iPO4=itot
      wrt_t  (itot)=.True.;       wrt_t_avg(itot)=.True.
      t_vname(itot)='PO4';        t_units  (itot)='mMol P m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Phosphate'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iNO3=itot
      wrt_t  (itot)=.True.;       wrt_t_avg(itot)=.True.
      t_vname(itot)='NO3';        t_units  (itot)='mMol N m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Nitrate'
      wrt_t_dia(itot)=.True.

      itot=itot+1;                iSIO3=itot
      wrt_t  (itot)=.False.;       wrt_t_avg(itot)=.True.
      t_vname(itot)='SiO3';       t_units  (itot)='mMol Si m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Silicate'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iNH4=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.True.
      t_vname(itot)='NH4';        t_units  (itot)='mMol N m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Ammonium'
      wrt_t_dia(itot)=.True.

      itot=itot+1;                iFE=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.True.
      t_vname(itot)='Fe';         t_units  (itot)='mMol Fe m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Iron'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iO2=itot
      wrt_t  (itot)=.True.;      wrt_t_avg(itot)=.True.
      t_vname(itot)='O2';         t_units  (itot)='mMol O2 m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Oxygen'
      wrt_t_dia(itot)=.True.

      itot=itot+1;                iDIC=itot
      wrt_t  (itot)=.True.;       wrt_t_avg(itot)=.True.
      t_vname(itot)='DIC';        t_units  (itot)='mMol C m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Dissolved inorganic carbon'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iALK=itot
      wrt_t  (itot)=.False.;       wrt_t_avg(itot)=.True.
      t_vname(itot)='ALK';        t_units  (itot)='mMol m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Alkalinity'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iDOC=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.True.
      t_vname(itot)='DOC';        t_units  (itot)='mMol C m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Dissolved organic carbon'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iDon=itot
      wrt_t  (itot) =.False.;      wrt_t_avg(itot)=.False.
      t_vname(itot)='DON';        t_units  (itot)='mMol N m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Dissolved organic nitrogen'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iDofe=itot
      wrt_t  (itot) =.False.;     wrt_t_avg(itot)=.False.
      t_vname(itot)='DOFE';       t_units  (itot)='mMol Fe m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Dissolved organic iron'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iDop=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.False.
      t_vname(itot)='DOP';        t_units  (itot)='mMol P m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Dissolved organic phosphorus'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iDopr=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.False.
      t_vname(itot)='DOPr';       t_units  (itot)='mMol P m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Refractory dissolved organic phosphorus'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iDonr=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.False.
      t_vname(itot)='DONr';       t_units  (itot)='mMol N m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Refractory dissolved organic nitrogen'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iZOOC=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.True.
      t_vname(itot)='zooC';       t_units  (itot)='mMol C m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Zooplankton'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iSPC=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.True.
      t_vname(itot)='spC';        t_units  (itot)='mMol C m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Small phytoplankton carbon'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iSPCHL=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.False.
      t_vname(itot)='spChl';      t_units  (itot)='mg Chl-a m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Small phytoplankton chlorophyll'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iSPFE=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.False.
      t_vname(itot)='spFe';       t_units  (itot)='mMol Fe m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Small phytoplankton iron'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iSPCACO3=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.False.
      t_vname(itot)='spCaCO3';    t_units  (itot)='mMol CaCO3 m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Small phytoplankton CaCO3'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iDIATC=itot
      wrt_t  (itot)=.False.;       wrt_t_avg(itot)=.True.
      t_vname(itot)='diatC';      t_units  (itot)='mMol C m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Diatom carbon'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iDIATCHL=itot
      wrt_t  (itot)=.False.;       wrt_t_avg(itot)=.False.
      t_vname(itot)='diatChl';    t_units  (itot)='mg Chl-a m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Diatom chlorophyll'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iDIATFE=itot
      wrt_t  (itot) =.False.;     wrt_t_avg(itot)=.False.
      t_vname(itot)='diatFe';     t_units  (itot)='mMol Fe m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Diatom Iron'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iDIATSI=itot
      wrt_t  (itot) =.False.;     wrt_t_avg(itot)=.False.
      t_vname(itot)='diatSi';     t_units  (itot)='mMol Si m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Diatom silicon'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iDiazc=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.True.
      t_vname(itot)='diazC';      t_units  (itot)='mMol C m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Diazotroph carbon'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iDiazchl=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.False.
      t_vname(itot)='diazChl';    t_units  (itot)='mg Chl-a m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Diazotroph chlorophyll'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iDiazfe=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.False.
      t_vname(itot)='diazFe';     t_units  (itot)='mMol Fe m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Diazotroph iron'
      wrt_t_dia(itot)=.False.

      !ntrc_bio_base=26 ! Total number of base bgc tracers. Hard-coded for now. Use itot later.

#ifdef Ncycle_SY
      itot=itot+1;                iNO2=itot
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.True.
      t_vname(itot)='NO2';        t_units  (itot)='mMol N m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Nitrite'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iN2=itot;
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.True.
      t_vname(itot)='N2';         t_units  (itot)='mMol N2 m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Dinitrogen'
      wrt_t_dia(itot)=.False.

      itot=itot+1;                iN2O=itot;
      wrt_t  (itot)=.False.;      wrt_t_avg(itot)=.True.
      t_vname(itot)='N2O';        t_units  (itot)='mMol N2O m-3'
      t_tname(itot)='';           t_ana_frc(itot)=1
      t_lname(itot)='Nitrous oxide'
      wrt_t_dia(itot)=.False.
#endif /* Ncycle_SY */

      ! total number of bgc tracers
