c*********************************************************************
c
c     program to extract .com files from LAMMPS simulations
c     for construct DME FORCEFIELD
c
c
c     author  A. Musetti July 2023
c
c*********************************************************************

      program mapeo
        parameter (mxatms=3000000)
        character*40 fname,linea
        character*500 commandline,AwkParte1,AwkParte2
        character*8 atmnam(mxatms)
        character*8 atmnamE(mxatms)
        character*8 aaa,bbb
c      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
        dimension x(mxatms),y(mxatms),z(mxatms)
        dimension IDS(mxatms),xS(mxatms),yS(mxatms),zS(mxatms)
        dimension IDE(mxatms),xE(mxatms),yE(mxatms),zE(mxatms)
        dimension IDA(mxatms)
        integer mxatmschk,count,NTOTS,nS,natms
        integer i,j,k,p,nC,nH,nO,nT,count2,NTOTc,kk,C2
        real(8) Rmax,xx,yy,zz,box,output

        NTOTc = 0
c      write(*,*)'Enter number of atoms of DME molecules in MD cell'
c      read(*,*) nS
        nS=725
c      write(*,*)'Maximum distance from Li to C (Rmax) ej. 15 A'
c      read(*,*) Rmax
        Rmax = 8
c      write(*,*)'Sphere around Li for DME (rmin) e.j. 2 A'
c      read(*,*) rmin

        natms = nS*16 + 2
c 1 Li y nS DME
        if(natms.gt.mxatms)then
          write(*,*)'Error - too many atoms in sphere radio. Max=',mxatms
          stop
        endif

c     CAMBIA LOS TIPOS POR NUMEROS
        CALL execute_command_line("sed -i -e 's/C /1 /g' xyz*.xyz" )
        CALL execute_command_line("sed -i -e 's/CT /2 /g' xyz*.xyz" )
        CALL execute_command_line("sed -i -e 's/H /3 /g' xyz*.xyz" )
        CALL execute_command_line("sed -i -e 's/O /4 /g' xyz*.xyz" )
        CALL execute_command_line("sed -i -e 's/Xx /5 /g' xyz*.xyz" )

c     Read xyz file
c      write(*,*)'Enter name of file xyz from traj LAMMPS'
c      read(*,*) fname
c        lee archivos xyz_30000.xyz   xyz_28000.xyz  xyz_26000.xyz ...

        do kk=1,5
        write(bbb,'(i0)') (1-kk) * 2000 + 30000
        fname = 'xyz_'//trim(bbb)//'.xyz'
        open(7,file=fname)
        AwkParte1 = "gawk '/^\s+"
        AwkParte2 = "/{a=$5} END{print a}' output > volumen"
        commandline = trim(AwkParte1)//trim(bbb)//trim(AwkParte2)
        call system(commandline)
        open(99,file='volumen')
        read(99, *) output
        box = output**(1./3.)
        close(99)
c        call system("rm volumen")
        read(7,*) mxatmschk
        if(mxatmschk.ne.natms)then
          write(*,*)'Error - check natms'
          write(*,*)'file ', mxatmschk,' natms', natms

          stop
        endif
        read(7,*) linea
c      NUEVO
c          read(7,*) linea
        do k=1,natms
        read(7,*) atmnam(k),x(k),y(k),z(k)
        enddo

        close(7)

c      write(*,*)'--------------------------------------------'
c       write(*,*)'                                            '
c       write(*,*)'           LEIDO ARCHIVO .xyz               '
c       write(*,*)'               natms =',natms
c       write(*,*)'                                            '
c       write(*,*)'--------------------------------------------'

c     Hasta aca ya lei todo el xyz
c     PRIMERO CENTRAMOS AL ULTIMO O

        do i=1,natms

c               Resta el ultimo O y pone en el centro
        xx = x(i) - x(natms) + box/2
        yy = y(i) - y(natms) + box/2
        zz = z(i) - z(natms) + box/2

        x(i) = xx
        y(i) = yy
        z(i) = zz

        if(xx.le.0)then
          x(i) = xx + box
        endif
        if(yy.le.0)then
          y(i) = yy + box
        endif
        if(zz.le.0)then
          z(i) = zz + box
        endif
        if(xx.gt.box)then
          x(i) = xx - box
        endif
        if(yy.gt.box)then
          y(i) = yy - box
        endif
        if(zz.gt.box)then
          z(i) = zz - box
        endif
        enddo


c       write(*,*)'--------------------------------------------'
c       write(*,*)'                                            '
c       write(*,*)'           CENTRADO O'
c       write(*,*)'           natoms =',natms,'                 '
c       write(*,*)'--------------------------------------------'

c      Ahora  a procesar
c     1ero leo y selecciono los C cercanos a la molec de O2
        count = 0
        countD = 0
        do i=1,natms
c           C:1
c            Selecciono solo el 1er C de los 2. CAMBIAR A +1
        nnn = countD*16+3
c       write(*,*)'nnn= ',nnn,' i= ',i
        if(i.eq.nnn)then
c               Resta el O ultimo
          xx = x(i) - x(natms)
          yy = y(i) - y(natms)
          zz = z(i) - z(natms)
          rS2 = xx*xx+yy*yy+zz*zz
          countD = countD + 1
c      write(*,*)'nnn= ',nnn,' rS2= ',rS2
          if(rS2.le.Rmax*Rmax)then
            count = count+1
            IDA(count) = i
            IDS(count) = count
            xS(count) = x(i)
            yS(count) = y(i)
            zS(count) = z(i)
c       write(*,*) count, IDA(count), IDS(count), xS(count), yS(count),
c     xzS(count)
          endif
        endif
        enddo
        NTOTS = count

c       write(*,*)'--------------------------------------------'
c       write(*,*)'                                            '
c       write(*,*)'           EXTRAIDOS C CERCANOS             '
c       write(*,*)'            TOTAL C = ', NTOTS,'            '
c       write(*,*)'                                            '
c       write(*,*)'--------------------------------------------'


        count = 1
        nO = 0
        nH = 0
        nC = 0
        nT = 0
c     Una vez extraidos los C_1, ahora busco todos los DME
c     o sea, los atomos alrededor

        do i=1,natms-2

        if(i.eq.IDA(count).and.atmnam(i).eq."1 ")then

          do j=1,16
          p = i+j
          count2 = (count-1)*16+j
          C2 = (count-1)*16+1
c          write(*,*) 'AAAAAAAAA', count2, C2, i
c             4:O
c                if(atmnam(p).eq."4 ")then
c                 write(*,*) 'para seguir',count,nC,nH,nO,nT
c                endif
c             4:O
          if(count2.eq.C2)then
c                if(atmnam(p).eq."4 ")then
            IDE(count2) = count
            atmnamE(count2)="C      "
            xE(count2) = x(p)
            yE(count2) = y(p)
            zE(count2) = z(p)
            nC=nC+1
          endif
c             2:CT
          if(atmnam(p).eq."2 ")then
            IDE(count2) = count
            atmnamE(count2)="CT     "
            xE(count2) = x(p)
            yE(count2) = y(p)
            zE(count2) = z(p)
            nT=nT+1
          endif
c              3:H
          if(atmnam(p).eq."3 ")then
            IDE(count2) = count
            atmnamE(count2)="H       "
            xE(count2) = x(p)
            yE(count2) = y(p)
            zE(count2) = z(p)
            nH=nH+1
          endif
c              4:O
          if(atmnam(p).eq."4 ")then
            IDE(count2) = count
            atmnamE(count2)="O       "
            xE(count2) = x(p-16)
            yE(count2) = y(p-16)
            zE(count2) = z(p-16)
            nO=nO+1
          endif
c             write(*,*) IDE(count2),atmnamE(count2),xE(count2),
c     xyE(count2), zE(count2), count, count2, j

          enddo
          count=count+1
        endif

        enddo

c       write(*,*)'--------------------------------------------'
c       write(*,*)'                                            '
c       write(*,*)'               EXTRAIDOS DME               '
c       write(*,*)'             TOTAL  = ', NTOTS,'            '
c       write(*,*)'            TOTAL C_1 = ', NTOTS,'            '
c       write(*,*)'            TOTAL C_2 = ', nC,'            '
c       write(*,*)'            TOTAL CT = ', nT,'            '
c       write(*,*)'            TOTAL H = ', nH,'            '
c       write(*,*)'            TOTAL O = ', nO,'            '
c       write(*,*)'                                            '
c       write(*,*)'--------------------------------------------'
c
c        if(nH.lt.3*nT)then
c       write(*,*)'--------------------------------------------'
c       write(*,*)'                                            '
c       write(*,*)'               ¡¡¡WARNING!!!!              '
c       write(*,*)'            TOTAL  = ', NTOTS,'            '
c       write(*,*)'           num H <   3 * num CT           '
c       write(*,*)'               AUMENTE R_min            '
c       write(*,*)'                                            '
c       write(*,*)'--------------------------------------------'
c        stop
c        endif
c        if(nH.gt.3*nT)then
c       write(*,*)'--------------------------------------------'
c       write(*,*)'                                            '
c       write(*,*)'               ¡¡¡WARNING!!!!              '
c       write(*,*)'            TOTAL  = ', NTOTS,'            '
c       write(*,*)'           num H >   3 * num C           '
c       write(*,*)'               DISMINUYA R_min            '
c       write(*,*)'                                            '
c       write(*,*)'--------------------------------------------'
c        stop
c        endif
c     Ahora a escribir el .com

      do j=1,NTOTS
        write(aaa,'(i0)') j + NTOTc
c       ! converting integer to string using a 'internal file'
             open(11,file='Archivo'//trim(aaa)//'.com')
             open(12,file='Archivo'//trim(aaa)//'.xyz')
             write(12,*) '18   '
             write(12,*) '    '
            write(11,'(a20)')'%Nproc=4            '
            write(11,'(a91)')"#n   b3lyp/6-311++g(2d,p) 5D Pop=None IOp(
     x3/124=30) SCF(maxcyc=3000) SCF=XQC CounterPoise=2"
            write(11,*) '    '
            write(11,*) '     '//trim(aaa)//' Preparing .coms file
     x  for dimer and fragment-pairs run'
            write(11,*) '    '
            write(11,*) '  0,3,  0,1,  0,3'
            write(11,'(a8,3f15.8,a2)') 'C       ',xS(j),yS(j),zS(j),' 1'
            write(12,'(a8,3f15.8)') 'C       ',xS(j),yS(j),zS(j)
            jj1 = (j-1) *  16  + 1
            jj2 = (j-1) *  16  + 15
c SIN H           jj1 = (j-1) *  3  + 1
c SIN H            jj2 = (j-1) *  3  + 3
              do k=jj1,jj2
            write(11,'(a2,3f15.8,a2)') atmnamE(k),xE(k),yE(k),zE(k),' 1'
            write(12,'(a2,3f15.8)') atmnamE(k),xE(k),yE(k),zE(k)
              enddo
          write(11,'(a8,3f15.8,a2)') 'O       ',x(natms-1),y(natms-1),
     xz(natms-1),' 2'
          write(12,'(a8,3f15.8)') 'O       ',x(natms-1),y(natms-1),
     xz(natms-1)
        write(11,'(a8,3f15.8,a2)') 'O       ',x(natms),y(natms),
     xz(natms),' 2'
          write(12,'(a8,3f15.8)') 'O       ',x(natms),y(natms),
     xz(natms)
           write(11,*) '    '
             close(11)
             close(12)
       enddo
        NTOTc = NTOTc + NTOTS
c      stop
        enddo
        return
        end

