	program coruij
************************************************************************
* Version 1.1 - 21 July 2011
*
* Usage:
*	coruij [-trim] [-noalt] [-main] [matrix] file1.pdb file2.pdb
*
*	-trim	Do not process atom pairs which differ by more than 1.5A
*	-noalt	Do not process atoms with an alternate conformation field
*	-alta	Treat conformation 'A' as being same as ' '
*	-strict Require matching atomname, chainid, res #, and res name
*		[default is to require only matching atom name and res #]
*	-iso	Second file is isotropic; report Biso, Beq, Best
*	-main	mainchain atoms only
*	[matrix]Rotation matrix to bring file2.pdb onto file1.pdb
*		several formats accepted:
*		-omat file.odb	O data block
*		-pmat file.xxx	PDB format rotation matrix
*
* In the current version every anisotropic atom in the first file
* must also be in the second file (in the same order).
*
* The PDB file convention for rotation/translation matrices is
* 	R11 R12 R13   T1
*	R21 R22 R23   T2
*	R31 R32 R33   T3
* that superimposes the atoms in file2.pdb onto those in file1.pdb.
*
* The O data block option has been tested using the output from lsqman
* A sample run goes something like this:
* lsqman
*	chain original
*	read m1 D.pdb
*	read m2 E.pdb
*	brute m1 D m2 E 50 25 50
*	show m1 m2
*	save m1 m2 DE.odb
*	quit
* coruij -trim -noalt -omat DE.odb D.pdb E.pdb 
*
* The program applies the transformation to the atomic coordinates of the
* second molecule, confirms that the molecules are superimposed 
* (rms diff in coords < 1.5A), then applies the same transformation to 
* the Uij of the second molecule.
*
* The program then writes to STDOUT a line for each atom giving the
* following quantities (as of 07-Jul-2011, subject to change!):
*   ATOM
*   RESTYP
*   RESNUM
*   coruij  (correlation coeff between un-normalized Uij, Vij
*   isocor1 (correlation between Uij and isotropic atom Ueq)
*   isocor2 (correlation between Vij and isotropic atom Veq)
*   corwij  (correlation between Uij and normalized Vij (force V'eq == Ueq))
*   SUV     corwij**2 / isocor1*isocor2
*   KLD     symmetrized Kullback-Leibler divergence
*   diffXYZ distance between paired atoms after superposition
*
*
************************************************************************
* To do:
*	- Should be forgiving of minor atom differences in the two files
*	- In particular should do better at ignoring alternate conformations
*	  in either file (there's code now, but I'm not sure it works right)
*	- Should optimize the transformation matrix before proceeding
*	- Is there a better normalization than multiplying by a constant
*	  to make Ueq = Veq?
*	- Is there a better "expected random value" than the correlation
*	  with an isotropic atom with Uiso = Ueq?
*	- should accept various formats for input matrix
*	  -pdbmat file.xxx   (as in current code)
*	  -omat   file.xxx   (O data block)
*         -xmat   file.xxx   (XtalView output format)
*	  [default]          assume coordinates already transformed
*
************************************************************************
* 24-Jun-1998	- Ethan A Merritt initial version
* 12-Nov-1998	- EAM add normalization of Uij matrices to make Ueq = Veq
* 15-Apr-1999	- EAM print out SUV
* 18-Apr-1999	-     fix bug in calculation of S(U,V)
* 15-May-2001	- print out column headings
* 07-Jul-2011	- Add Kullback-Leibler divergence
* 19-Jul-2011	- -iso option to evaluate Biso / Beq / Best
*
************************************************************************
	implicit none

c external routines
	integer  iargc
	external inv3x3, trn3x3
	real     inv3x3

c I/O units
	integer  in1, in2, in3, stdout, stderr
	parameter (in1=1, in2=2, in3=3, stdout=6, stderr=0)
	character*80 card
	character*64 filename1, filename2, filename3
	character*2  main

c Command line options
	character*64 flag
	logical      trimflag, altflag, mainflag, altaflag
	logical      ahead2, strictflag, isoflag

c atom arrays
	integer    MAXATOM
	parameter (MAXATOM=10000)
	character*28 atom(MAXATOM)
	real         x1(MAXATOM), y1(MAXATOM), z1(MAXATOM)
	real         x2(MAXATOM), y2(MAXATOM), z2(MAXATOM)
	real         xt(MAXATOM), yt(MAXATOM), zt(MAXATOM)
	real         Uij_1(6,MAXATOM)
	real         Uij_2(6,MAXATOM)
	real         diffxyz(MAXATOM)
	real         uijcor(MAXATOM), isocor1, isocor2
	real         corwij, Suv
	real         Ueq, Veq

c atom specifier components that we use to match atoms
	character*4  atname1,  atname2
	character*3  resname1, resname2
	character*1  altloc1,  altloc2
	character*1  chainid1, chainid2
	integer      resseq1,  resseq2

c transformation matrices
	real TRAN(3)
	real TMAT(3,3), TMATT(3,3)
	real TINV(3,3), TINVT(3,3)

c Uij matrices and determinants
	real U(3,3),  UINV(3,3)
	real V(3,3),  VINV(3,3)
	real W(3,3),  WINV(3,3)
	real          ZINV(3,3)
	real VV(3,3), VVINV(3,3)
	real KLD, KL1(3,3), KL2(3,3), KL3(3,3)
	real det, det_U, det_V, det_W, det_VV
	real rescale

c Keep statistics for spatial fit
	real del, sumdel, sumdel2, rmsdel, avgdel

c Variables used by the -iso comparison mode
	real Beq, Biso, Best, eigen(3), Anisotropy

c General bookkeeping
	integer i, j, k, iatm
	integer natom, nskip2
	integer ntoofar
	integer njunk

c say hello
	write(stderr,'(A58)') 
     &     '    *****************************************************',
     &     '    *  coruij V1.1                          21-Jul-2011 *',
     &     '    *                                                   *',
     &     '    *  comments and suggestions to:                     *',
     &     '    *      Ethan A Merritt  -  merritt@u.washington.edu *',
     &     '    *****************************************************'

c Initialize
	sumdel = 0
	sumdel2 = 0

c Assume identity transformation
	do i = 1,4
	do j = 1,4
	   TMAT(i,j) = 0.0
	enddo
	TMAT(i,i) = 1.0
	enddo
	do i = 1,3
	   TRAN(i) = 0.0
	enddo

c Attempt to open input files
	if (iargc().lt.2) goto 900
	trimflag = .false.
	altflag  = .false.
	altaflag = .false.
	mainflag = .false.
	ahead2   = .false.
	strictflag = .false.
	isoflag  = .false.
	i = 0

  10	continue
	call getarg(i+1,flag)

	if (flag(1:5).eq.'-help') goto 900
	
	if (flag(1:5).eq.'-trim') then
	    trimflag = .true.
	    i = i + 1
	    goto 10
	endif

	if (flag(1:6).eq.'-alta') then
	    altaflag = .true.
	    i = i + 1
	    goto 10
	endif

	if (flag(1:6).eq.'-noalt') then
	    altflag = .true.
	    i = i + 1
	    goto 10
	endif

	if (flag(1:5).eq.'-main') then
	    mainflag = .true.
	    i = i + 1
	    goto 10
	endif

	if (flag(1:7).eq.'-strict') then
	    strictflag = .true.
	    i = i + 1
	    goto 10
	endif

	if (flag(1:4).eq.'-iso') then
	    isoflag = .true.
	    i = i + 1
	    goto 10
	endif

	if (iargc().lt.i+2) goto 900
	if (flag(1:5).eq.'-pmat') then
	    call getarg( i+2, filename3)
	    open ( unit=in3, file=filename3, status='old', err=903 )
	    read ( in3, *, err=904 ) (TMAT(1,j),j=1,3),TRAN(1)
	    read ( in3, *, err=904 ) (TMAT(2,j),j=1,3),TRAN(2)
	    read ( in3, *, err=904 ) (TMAT(3,j),j=1,3),TRAN(3)
	    i = i + 2
	    goto 10
	endif
	if (flag(1:5).eq.'-omat') then
	    call getarg( i+2, filename3)
	    open ( unit=in3, file=filename3, status='old', err=903 )
	    read ( in3, '(A80)', err=904 ) card
	    write (stderr,*) card
	    read ( in3, '(A80)', err=904 ) card
	    write (stderr,*) card
	    read ( in3, *, err=904 ) (TMAT(j,1),j=1,3)
	    read ( in3, *, err=904 ) (TMAT(j,2),j=1,3)
	    read ( in3, *, err=904 ) (TMAT(j,3),j=1,3)
	    read ( in3, *, err=904 ) (TRAN(j),j=1,3)
	    i = i + 2
	    goto 10
	endif

	call getarg( i+1, filename1 )
	open ( unit=in1, file=filename1, status='old', err=901 )
	call getarg( i+2, filename2 )
	open ( unit=in2, file=filename2, status='old', err=902 )
	goto 100

c Error reporting
900	write(stderr,*) 'Usage: ',
     &    'coruij [-trim] [-noalt] [-strict] [-main] [matrix] ',
     &    'file1.pdb file2.pdb'
        write(stderr,*) '              [matrix] = -omat file.odb'
        write(stderr,*) '                      or -pmat file.matrix'
     	write(stderr,*) ' output columns (subject to change!):'
     	write(stderr,'((5x,A))')
     & 	' 1  ATOM',
     & 	' 2  RESTYP',
     & 	' 3  RESNUM',
     & 	' 4  coruij  (correlation coeff between un-normalized Uij, Vij',
     & 	' 5  isocor1 (correlation between Uij and isotropic atom Ueq)',
     & 	' 6  isocor2 (correlation between Vij and isotropic atom Veq)',
     & 	' 7  corwij  (correlation between Uij and normalized Vij)',
     & 	' 8  SUV     = corwij**2 / isocor1*isocor2',
     & 	' 9  KLD     symmetrized Kullback-Leibler divergence',
     &  '10  diffXYZ distance between paired atoms after superposition'
     	write(stderr,*) ' In -iso mode the output is:'
	write(stderr,'((5x,A))')
     &  '    ATOM RESTYP RESNUM Anisotropy Biso Beq Best'
	call exit(-1)
901	write(stderr,*) 'Cannot open ',filename1
	call exit(-1)
902	write(stderr,*) 'Cannot open ',filename2
	call exit(-1)
903	write(stderr,*) 'Cannot open ',filename3
	call exit(-1)
904	write(stderr,*) 'Cannot read TMAT from ',filename3
	call exit(-1)
905	write(stderr,*) 'problem reading from ',filename1
	write(stderr,*) 'current line is: ',card
	call exit(-1)
906	write(stderr,*) 'problem reading from ',filename2
	write(stderr,*) 'current line is: ',card
	call exit(-1)

c
  100	continue
  	write (stderr,101) (TMAT(1,i),i=1,3),TRAN(1),
     &                     (TMAT(2,i),i=1,3),TRAN(2),
     &                     (TMAT(3,i),i=1,3),TRAN(3)
  101	format(' Input transformation matrix:',3(/,3f9.4,f12.4),/)

c Walk through the two PDB input files in parallel
c First we find an atom in file1 with ANISOU records,
c then we look for a corresponding atom in file2
c
	natom   = 0
	nskip2  = 0
	ntoofar = 0
c
  200	continue
  	read (in1,'(A80)',end=300) card
	if (card(1:4).ne.'ATOM' .and. card(1:6).ne.'HETATM') goto 200
	if (altflag .and. card(17:17).ne.' ') goto 200
	if (altaflag .and. card(17:17).ne.' ' .and. card(17:17).ne.'A') 
     &	    goto 200
	main = card(14:15)
	if (mainflag .and. .not.
     &	   (main.eq.'C ' .or. main.eq.'N ' .or. main.eq.'O ' .or.
     &      main.eq.'CA')) goto 200
	natom = natom + 1
	atom(natom)(1:27) = card(1:27)
	read (card,201,err=905) x1(natom),y1(natom),z1(natom)
  201	format(30x,3f8.3)
c
	atname1  = card(13:16)
	resname1 = card(18:20)
	altloc1  = card(17:17)
	chainid1 = card(22:22)
	read (card(23:26),'(I4)') resseq1
c
	read (in1,'(A80)') card
	if (card(1:6).ne.'ANISOU') then
	    natom = natom - 1
	    goto 200
	endif
	read (card(29:70),202,err=905) (Uij_1(i,natom),i=1,6)
  202	format(6f7.0)
  	do i=1,6
	    Uij_1(i,natom) = Uij_1(i,natom) * 0.0001
	enddo

c Enforce match between atoms in the two files
c (default)	atom name and residue # must match; chain ID and resnam ignored
c  -strict	all entries must match
c
  210	continue
  	if (ahead2) goto 212
	read (in2,'(A80)',end=299) card
	if (card(1:4).ne.'ATOM' .and. card(1:6).ne.'HETATM') goto 210
	if (altflag .and. card(17:17).ne.' ') goto 210
	if (altaflag .and. card(17:17).ne.' ' .and. card(17:17).ne.'A') 
     &	    goto 210
	main = card(14:15)
	if (mainflag .and. .not.
     &	   (main.eq.'C ' .or. main.eq.'N ' .or. main.eq.'O ' .or.
     &      main.eq.'CA')) goto 210
c
	atname2  = card(13:16)
	resname2 = card(18:20)
	altloc2  = card(17:17)
	chainid2 = card(22:22)
	read (card(23:26),'(I4)') resseq2
c
	read (card,201,err=906) x2(natom),y2(natom),z2(natom)
	if (isoflag) then
	    read (card(61:66),'(F6.0)') Uij_2(1,natom)
	else
	    read (in2,'(A80)') card
	    if (card(1:6).eq.'ANISOU') then
		read (card(29:70),202,err=906) (Uij_2(i,natom),i=1,6)
		do i=1,6
		    Uij_2(i,natom) = Uij_2(i,natom) * 0.0001
		enddo
	    else
		atname2 = 'XXXX'
	    endif
	endif
  212	continue
  	ahead2 = .false.
c
c	Check if chain in file 1 is different length than same chain in file 2
	if (chainid1.ne.chainid2 .and. strictflag) then
	    if (resseq1.lt.resseq2) then
	    	nskip2 = nskip2 + 1
		write (stderr,*) '   skipping ',
     &				atname2,resname2,altloc2,chainid2,resseq2
	    	goto 210
	    endif
	    if (resseq1.gt.resseq2) then
		ahead2 = .true.
	    	natom = natom - 1
		goto 200
	    endif
	endif
c
c	Check if residues in the two files have gotten out of sync
	if (resseq1.lt.resseq2) then
	    ahead2 = .true.
	    natom = natom - 1
	    goto 200
	endif
	if (resseq1.gt.resseq2) then
	    nskip2 = nskip2 + 1
	    write (stderr,*) '   skipping ',
     &			atname2,resname2,altloc2,chainid2,resseq2
	    goto 210
	endif
c
c	Check if the atom names match
	if (atname1.ne.atname2) then
	    nskip2 = nskip2 + 1
	    write (stderr,*) '   skipping ',
     &			atname2,resname2,altloc2,chainid2,resseq2
	    goto 210
	endif
c
c	Check if the alternate conformation flags match
	if ( altloc1.eq.' ' .and. altloc2.eq.'A' 
     &  .or. altloc1.eq.'A' .and. altloc2.eq.' ' ) then
     	    njunk = njunk + 1
c
        else if (altloc1.ne.altloc2) then
	    nskip2 = nskip2 + 1
	    write (stderr,*) '   skipping ',
     &			atname2,resname2,altloc2,chainid2,resseq2
	    goto 210
	endif
c
c	Only check residue name if -strict flag is given
	if (strictflag) then
	    if (resname1.ne.resname2) then
	    	nskip2 = nskip2 + 1
		write (stderr,*) '   skipping ',
     &				atname2,resname2,altloc2,chainid2,resseq2
	    	goto 210
	    endif
	endif
c
c	Accept match
	goto 200
c
  299	natom = natom - 1

c All done reading atoms
  300	continue
  	write (stderr,301) natom
  301	format(' Found',i8,' matching atoms')
	if (nskip2.gt.0) write (stderr,302) nskip2
  302	format(' Had to skip',i4,' atoms in 2nd file that did not match')

c
c	  Cosmetic changes to atom identifier for the sake of sorting
c	  We will force there to be exactly three entities printed.
c	  PDB format is just a mess:
c	  cols 13:16	atom
c	  col     17	alternate conf
c	  cols 18:20	residue
c	  col     22	chain
c	  cols 23:27	resnum
c
	do iatm = 1, natom

	    do i = 16, 13, -1
		if (ATOM(iatm)(i:i) .ne. ' ') j = i
	    enddo
	    do i = 13, 17
		if (ATOM(iatm)(i:i) .ne. ' ') k = i
	    enddo
	    do i = j, k
		if (ATOM(iatm)(i:i) .eq. ' ') ATOM(iatm)(i:i) = '_'
	    enddo
c	    if (ATOM(iatm)(17:17) .ne. ' ') then
c	      do i = 18,19
c		if (ATOM(iatm)(i:i) .eq. ' ') ATOM(iatm)(i:i) = '_'
c	      enddo
c	    endif
	    if (ATOM(iatm)(22:22) .ne. ' ') then
	      do i = 23,25
		if (ATOM(iatm)(i:i) .eq. ' ') ATOM(iatm)(i:i) = '_'
	      enddo
	    endif
	enddo


c Apply TMAT to 2nd set of coordinates
	do i = 1, natom
	    xt(i) = x2(i)*TMAT(1,1) + y2(i)*TMAT(1,2) + z2(i)*TMAT(1,3)
     &            + TRAN(1)
	    yt(i) = x2(i)*TMAT(2,1) + y2(i)*TMAT(2,2) + z2(i)*TMAT(2,3)
     &            + TRAN(2)
	    zt(i) = x2(i)*TMAT(3,1) + y2(i)*TMAT(3,2) + z2(i)*TMAT(3,3)
     &            + TRAN(3)
	enddo
	do i = 1, natom
	    del = (x1(i)-xt(i))**2 + (y1(i)-yt(i))**2 + (z1(i)-zt(i))**2
	    diffxyz(i) = sqrt(del)
	    sumdel2 = sumdel2 + del
	    sumdel  = sumdel  + diffxyz(i)
	    if (diffxyz(i) .gt. 1.5) ntoofar = ntoofar + 1
	enddo
	avgdel = sumdel / float(natom)
	rmsdel = sqrt( sumdel2 / float(natom) )
	write (stderr,*) 'RMS difference in coordinates after ',
     &                   'superposition: ', rmsdel
     	if (rmsdel .gt. 1.5) then
	    write (stderr,*) 
     &            '>>> That seems a bit high. Are you sure ',
     &            'the matrix is right? <<<'
	else
	    write (stderr,303) ntoofar
  303	    format(i8,' atom pairs farther apart than cutoff of 1.5A')
	    if (trimflag) write (stderr,304)
  304	    format(8x,' and will be skipped (-trim option)')
	endif
	    
c Should optimize TMAT and transformation here
c

c Construct various transposes and inverses of rotation matrix
c
	det = inv3x3( TINV,  TMAT )
	call assert( abs(det-1.0).lt.0.01, 'TMAT determinant != 1' )
	call trn3x3( TMATT, TMAT )
	call trn3x3( TINVT, TINV )

c OK, here we go.
      if (.not.isoflag) goto 310
c The isotropic comparison mode
c Loop over all the atoms and track agreement in Biso / Beq / Best
c
      do i = 1, natom

c Unpack Uij for the aniso version of this atom
c
	U(1,1) = Uij_1(1,i)
	U(2,2) = Uij_1(2,i)
	U(3,3) = Uij_1(3,i)
	U(1,2) = Uij_1(4,i)
	U(2,1) = Uij_1(4,i)
	U(1,3) = Uij_1(5,i)
	U(3,1) = Uij_1(5,i)
	U(2,3) = Uij_1(6,i)
	U(3,2) = Uij_1(6,i)

	det_U = inv3x3( Uinv, U )
	if (det_U .le. 0) then
	    write (stderr,333) atom(i)(13:17),atom(i)(18:27)
	    goto 309
	endif

c Get the Eigenvalues for it
c
	call eigen3x3(U,eigen)
	Best = sqrt( (eigen(1) + eigen(2) + eigen(3))
     &            / (1/eigen(1) + 1/eigen(2) + 1/eigen(3)))
	Best = Best * 8. * 3.14159**2
	Anisotropy = eigen(3) / eigen(1)
	Beq = (U(1,1) + U(2,2) + U(3,3)) / 3.0
	Beq = Beq * 8. * 3.14159**2
	Biso = Uij_2(1,i)

c Write entry for this pair into output file
c
	write (stdout,305) atom(i)(13:17),atom(i)(18:27),
     &            Anisotropy,
     &            Biso, Beq, Best
  305	  format(A5,1X,A10,F10.4,4X,3F10.4,4X,3F10.4)
  309 continue

      enddo

      goto 400


  310 continue
c The default (anisotropic) comparison mode
c Loop over all the atoms and track agreement in Uij's
c
      do i = 1, natom

c Unpack Uij of 1st set of atoms
c
	U(1,1) = Uij_1(1,i)
	U(2,2) = Uij_1(2,i)
	U(3,3) = Uij_1(3,i)
	U(1,2) = Uij_1(4,i)
	U(2,1) = Uij_1(4,i)
	U(1,3) = Uij_1(5,i)
	U(3,1) = Uij_1(5,i)
	U(2,3) = Uij_1(6,i)
	U(3,2) = Uij_1(6,i)

c Apply rotation matrix to Uij of 2nd set of atoms
c V' = TINVT * V * TINV
c
	V(1,1) = Uij_2(1,i)
	V(2,2) = Uij_2(2,i)
	V(3,3) = Uij_2(3,i)
	V(1,2) = Uij_2(4,i)
	V(2,1) = Uij_2(4,i)
	V(1,3) = Uij_2(5,i)
	V(3,1) = Uij_2(5,i)
	V(2,3) = Uij_2(6,i)
	V(3,2) = Uij_2(6,i)
	call mul3x3( W, V, TINV )
	call mul3x3( V, TINVT, W )

c Take inverses of the two matrices
c
	det_U = inv3x3( Uinv, U )
	det_V = inv3x3( Vinv, V )

	if (det_U .le. 0 .or. det_V .le. 0) then
	    write (stderr,333) atom(i)(13:17),atom(i)(18:27)
	    goto 399
	endif
  333	format(' Non-positive definite matrix for ',a5,1x,a10)

c And their sum (for the covariance)
c
	call add3x3( Winv, Uinv, Vinv )
	det_W = 1. / inv3x3( W, Winv )

c Now for the correlation coefficient
c
	uijcor(i) = (8.0 * det_W) / sqrt( det_U * det_V )
	uijcor(i) = sqrt(uijcor(i))

c Check distance separating this atom pair
c
	if (trimflag .and. diffxyz(i).gt.1.5) goto 399

c Get some idea of baseline correlation by comparing to
c isotropic atom with same Beq
c
	Ueq = (U(1,1)+U(2,2)+U(3,3)) / 3.0
	do j = 1, 3
	    Zinv(j,j) = 1./Ueq
	enddo
	call add3x3( Winv, Uinv, Zinv )
	det_W = 1. / inv3x3( W, Winv )
	isocor1 = (8.0 * det_W) / sqrt( Ueq*Ueq*Ueq * det_U )
	isocor1 = sqrt(isocor1)
c
	Veq = (V(1,1)+V(2,2)+V(3,3)) / 3.0
	do j = 1, 3
	    Zinv(j,j) = 1./Veq
	enddo
	call add3x3( Winv, Vinv, Zinv )
	det_W = 1. / inv3x3( W, Winv )
	isocor2 = (8.0 * det_W) / sqrt( Veq*Veq*Veq * det_V )
	isocor2 = sqrt(isocor2)

c 10-Nov-1998
c Re-scale 2nd atom so that it has the same Beq as the 1st, 
c but keeps its anisotropy
c
	rescale = (Ueq/Veq)
	do j = 1, 3
	do k = 1, 3
	    VV(j,k) = V(j,k) * rescale
	enddo
	enddo
	det_VV = inv3x3( VVinv, VV )
	
	call add3x3( Winv, Uinv, VVinv )
	det_W  = 1. / inv3x3( W, Winv )
	corwij = (8.0 * det_W) / sqrt( det_U * det_VV )
	corwij = sqrt(corwij)


c 15-Apr-1999
c Define Suv = corwij**2 / isocor1*isocor2
	Suv   = (corwij*corwij) / (isocor1*isocor2)

c 07-Jul-2011
c And the symmetrized Kullback-Leibler divergence
C This formula is taken from Garib Murshudov's teaching notes
	call mul3x3(KL1, U, Vinv)
	call mul3x3(KL2, V, Uinv)
	call add3x3(KL3, KL1, KL2)
	KLD = KL3(1,1) + KL3(2,2) + KL3(3,3) - 6.0

c Write entry for this pair into output file
c
	write (stdout,306) atom(i)(13:17),atom(i)(18:27),
     &            uijcor(i), isocor1, isocor2,
     &            corwij, Suv
     &            , KLD
     &            , diffxyz(i)
  306	  format(A5,1X,A10,8F10.4)

c That's it for the grand loop over atoms
c
  399	continue
      enddo

c Now a statistical summary at the end
c
  400	continue

c
	end



c***********************************************************************
c Here follow some matrix manipulation routines and the like           *
c that probably should live in a library somewhere                     *
c***********************************************************************

	subroutine trn3x3( A, B )
	real    A(3,3), B(3,3)
	integer i, j
	do i = 1, 3
	do j = 1, 3
	    A(i,j) = B(j,i)
	enddo
	enddo
	return
	end

	subroutine mul3x3( C, A, B )
	real A(3,3), B(3,3), C(3,3)
	integer i,j
	do i = 1,3
	do j = 1,3
	    C(i,j) = A(i,1)*B(1,j) + A(i,2)*B(2,j)
     &		   + A(i,3)*B(3,j)
	enddo
	enddo
	return
	end

	subroutine add3x3( C, A, B )
	real A(3,3), B(3,3), C(3,3)
	integer i,j
	do i = 1,3
	do j = 1,3
	    C(i,j) = A(i,j) + B(i,j)
	enddo
	enddo
	return
	end
	
	function inv3x3( A, B )
	real    inv3x3
	real    A(3,3), B(3,3)
	real    tmp(3,3), D
	integer index(3)
	integer i, j
	do i=1,3
	    do j=1,3
		TMP(i,j) = B(i,j)
		A(i,j) = 0.
	    enddo
	    A(i,i) = 1.
	enddo
	call ludcmp( TMP, 3, index, D )
	inv3x3 = D * TMP(1,1)*TMP(2,2)*TMP(3,3)
	do j=1,3
	    call lubksb( TMP, 3, index, A(1,j) )
	enddo
	return
	end

	subroutine assert (logic, message)
	logical logic
	character*(*) message
	if (logic) return
	write(0,*) '>>> ',message,' <<<'
	stop
	end


************************************************************************
*              Matrix inversion via LU decomposition                   *
*	adapted from Numerical Recipes in Fortran (1986)               *
************************************************************************
*
CCC	input  NxN matrix A is replaced by its LU decomposition
CC	output index(N) records row permutation due to pivoting
C	output D is parity of row permutations 
c
	subroutine ludcmp( A, n, index, D )
	implicit  NONE
	integer   n,index(n)
	real      A(n,n)
	real      D
c
	integer   i,imax,j,k
	real      aamax,dum,sum
	integer    NMAX
	parameter (NMAX=10)
	real      vv(NMAX)
c
	d = 1.
	do i=1,n
	    aamax = 0.
	    do j=1,n
	    	if (abs(A(i,j)).gt.aamax) aamax = abs(A(i,j))
	    enddo
	    call assert(aamax.ne.0.,'Singular matrix')
	    vv(i) = 1. / aamax
	enddo
c
	do j=1,n
	    if (j.gt.1) then
		do i=1,j-1
		    sum = A(i,j)
		    if (i.gt.1) then
			do k=1,i-1
			    sum = sum - A(i,k)*A(k,j)
			enddo
			A(i,j) = sum
		    endif
		enddo
	    endif
	    aamax = 0.
	    do i=j,n
		sum = A(i,j)
		if (j.gt.1) then
		    do k=1,j-1
			sum = sum - A(i,k)*A(k,j)
		    enddo
		    A(i,j) = sum
		endif
		dum = vv(i) * abs(sum)
		if (dum.ge.aamax) then
		    imax = i
		    aamax = dum
		endif
	    enddo
	    if (j.ne.imax) then
		do k=1,n
		    dum = A(imax,k)
		    A(imax,k) = A(j,k)
		    A(j,k) = dum
		enddo
		d = -d
		vv(imax) = vv(j)
	    endif
	    index(j) = imax
	    call assert(A(j,j).ne.0.,'Singular matrix')
	    if (j.ne.n) then
		dum = 1. / A(j,j)
		do i=j+1,n
		    A(i,j) = A(i,j) * dum
		enddo
	    endif
	enddo
	return
	end

CCC	corresponding back-substitution routine
CC
C
	subroutine lubksb( A, N, index, B )
	implicit NONE
	integer  n, index(n)
	real     A(n,n), B(n)
c
	integer  i,ii,j,ll
	real     sum
c
	ii = 0
	do i=1,n
	    ll = index(i)
	    sum = B(ll)
	    B(ll) = B(i)
	    if (ii.ne.0) then
		do j=ii,i-1
		    sum = sum - A(i,j)*B(j)
		enddo
	    else if (sum.ne.0.) then
		ii = i
	    endif
	    B(i) = sum
	enddo
c
	do i=n,1,-1
	    sum = B(i)
	    if (i.lt.n) then
		do j=i+1,n
		    sum = sum - A(i,j)*B(j)
		enddo
	    endif
	    B(i) = sum / A(i,i)
	enddo
	return
	end


************************************************************************
* Find eigenvalues of a 3x3 symmetric matrix M                         *
* Non-destructive for M.  Sorted eigenvalues returned in E.            *
* Reference: Oliver K. Smith (1961) CACM 4:168                         *
* Note: this method does not work for singular matrices.               *
************************************************************************
CCC
CC
C
	subroutine eigen3x3( M, E )
	real M(3,3)
	real E(3)
C
	real K(3,3)
	real trace, det, p, q, phi

	trace = (M(1,1) + M(2,2) + M(3,3)) / 3.
	do i = 1,3
	  do j = 1,3
	    K(i,j) = M(i,j)
	  enddo
	  K(i,i) = K(i,i) - trace
	enddo

	det = K(1,1)*K(2,2)*K(3,3)
     c      + K(1,2)*K(2,3)*K(3,1) + K(1,3)*K(2,1)*K(3,2)	
     c      - K(1,1)*K(2,3)*K(3,2) - K(3,3)*K(1,2)*K(2,1)
     c      - K(1,3)*K(2,2)*K(3,1)

	q = det / 2.
	p = 0
	do i = 1,3
	  do j = 1,3
	    p = p + K(i,j)*K(i,j)
	  enddo
	enddo
	p = p / 6.

	phi = acos(q / p**(1.5)) / 3.
 
C	Check for round-off error
	if (abs(q) .ge. abs(p**(1.5))) phi = 0.

	E(1) = trace + 2. * sqrt(p) * cos(phi)
	E(2) = trace - sqrt(p) * (cos(phi) + sqrt(3.)*sin(phi))
	E(3) = trace - sqrt(p) * (cos(phi) - sqrt(3.)*sin(phi))

	if (E(2).gt.E(1)) then
	  p = E(1)
	  E(1) = E(2)
	  E(2) = p
	endif
	if (E(3).gt.E(2)) then
	  p = E(2)
	  E(2) = E(3)
	  E(3) = p
	endif

	return
	end
