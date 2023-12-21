MODULE gc_mod
    USE core_mod
    USE simdfunctions_mod, ONLY: divide0
    USE ibmodel_mod, ONLY: ibmodel_t
    USE ibconst_mod, ONLY: maccur
    USE noib_mod, ONLY: noib_t
    USE bubvbw_mod, ONLY: bubvbw
    USE gc_restrict_mod, ONLY: gc_restrict_t
    USE gc_blockbp_mod, ONLY: gc_blockbp_t
    USE gc_stencils_mod, ONLY: gc_stencils_t

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE, EXTENDS(noib_t) :: gc_t
        INTEGER(intk) :: ncells
        INTEGER(intk), ALLOCATABLE :: icells(:)
        INTEGER(intk), ALLOCATABLE :: icellspointer(:)

        REAL(realk), ALLOCATABLE :: xpsw(:, :), nvecs(:, :), ucell(:, :)
        REAL(realk), ALLOCATABLE :: arealist(:)
        INTEGER(intk), ALLOCATABLE :: bodyid(:)
        INTEGER(intk), ALLOCATABLE :: bzelltyp(:)

        TYPE(gc_stencils_t) :: stencils
        CHARACTER(len=mglet_filename_max) :: stencilfile
        CHARACTER(len=mglet_filename_max), ALLOCATABLE :: stlnames(:)
    CONTAINS
        PROCEDURE :: blockbp
        PROCEDURE :: read_stencils
        PROCEDURE :: divcal
        PROCEDURE, PRIVATE :: calc_nvecs
        PROCEDURE, NOPASS, PRIVATE :: calc_nvecs_grid
        FINAL :: destructor
    END TYPE gc_t

    INTERFACE gc_t
        MODULE PROCEDURE :: constructor
    END INTERFACE gc_t

    PUBLIC :: gc_t, constructor

CONTAINS
    FUNCTION constructor() RESULT(ib)
        ! Subroutine arguments
        CLASS(ibmodel_t), ALLOCATABLE :: ib

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: bp(:), bu(:), bv(:), bw(:)
        REAL(realk), POINTER, CONTIGUOUS :: au(:), av(:), aw(:)

        IF (myid == 0) THEN
            WRITE(*, '(A)') "Using 'ghostcell' immersed boundary method"
            WRITE(*, '()')
        END IF

        ALLOCATE(gc_t :: ib)
        ALLOCATE(gc_restrict_t :: ib%restrict_op)
        ALLOCATE(gc_blockbp_t :: ib%blockbp_op)

        ib%type = "GHOSTCELL"

        CALL set_field("BP", dwrite=.TRUE.)
        CALL set_field("BU")
        CALL set_field("BV")
        CALL set_field("BW")

        CALL set_field("AU")
        CALL set_field("AV")
        CALL set_field("AW")

        CALL set_field("SDIV")

        CALL get_fieldptr(bp, "BP")
        CALL get_fieldptr(bu, "BU")
        CALL get_fieldptr(bv, "BV")
        CALL get_fieldptr(bw, "BW")

        CALL get_fieldptr(au, "AU")
        CALL get_fieldptr(av, "AV")
        CALL get_fieldptr(aw, "AW")

        bp = 1.0
        bu = 1.0
        bv = 1.0
        bw = 1.0

        au = 1.0
        av = 1.0
        aw = 1.0

        SELECT TYPE(ib)
        TYPE IS (gc_t)
            CALL fort7%get_value("/ib/stencilfile", ib%stencils%file, &
                "ib_stencils.h5")

            ALLOCATE(ib%icells(ngrid))
            ALLOCATE(ib%icellspointer(ngrid))
            ib%icells = 0
            ib%icellspointer = 0
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END FUNCTION constructor


    SUBROUTINE destructor(this)
        TYPE(gc_t), INTENT(inout) :: this
        IF (ALLOCATED(this%restrict_op)) THEN
            DEALLOCATE(this%restrict_op)
        END IF
    END SUBROUTINE destructor


    SUBROUTINE blockbp(this)
        ! Subroutine arguments
        CLASS(gc_t), INTENT(inout) :: this

        SELECT TYPE(blockbp_op => this%blockbp_op)
        TYPE IS (gc_blockbp_t)
            CALL blockbp_op%blockbp(this%icells, this%icellspointer, &
                this%stencils)
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE blockbp


    SUBROUTINE read_stencils(this)
        ! Subroutine arguments
        CLASS(gc_t), INTENT(inout) :: this

        ! Local variables
        TYPE(field_t), POINTER :: bp, bu, bv, bw
        TYPE(field_t), POINTER :: au, av, aw

        CALL this%stencils%read()

        CALL get_field(bp, "BP")
        CALL this%stencils%get_bp(bp)

        CALL get_field(bu, "BU")
        CALL get_field(bv, "BV")
        CALL get_field(bw, "BW")
        CALL bubvbw(bp, bu, bv, bw)

        CALL get_field(au, "AU")
        CALL get_field(av, "AV")
        CALL get_field(aw, "AW")
        CALL this%stencils%get_auavaw(bu, bv, bw, au, av, aw)

        CALL this%stencils%get_icells(this%icells)
        CALL this%blockbp_op%set_icellspointer(this%icells, this%icellspointer)

        this%ncells = SUM(this%icells)
        ALLOCATE(this%xpsw(3, this%ncells))
        ALLOCATE(this%ucell(3, this%ncells))
        ALLOCATE(this%bodyid(this%ncells))
        ALLOCATE(this%arealist(this%ncells))
        ALLOCATE(this%bzelltyp(idim3d))
        CALL this%stencils%get_intersected(this%icells, this%icellspointer, &
            this%bodyid, this%xpsw, this%ucell, this%bzelltyp)

        ALLOCATE(this%nvecs(4, this%ncells))
        CALL this%calc_nvecs(this%bzelltyp, au, av, aw, this%icells, &
            this%icellspointer, this%xpsw, this%nvecs, this%arealist)

        ALLOCATE(this%stlnames, SOURCE=this%stencils%stlnames)

        CALL this%stencils%finish()
    END SUBROUTINE read_stencils

    SUBROUTINE giteig(this)
        ! Subroutine arguments
        CLASS(gc_t), INTENT(inout) :: this

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: igr, igrid
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        TYPE(field_t), POINTER :: au_f, av_f, aw_f

        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: gsae(:,:,:),  gsaw(:,:,:), &
            gsan(:,:,:), gsas(:,:,:), gsat(:,:,:), gsab(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: gsap(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: au(:, :, :), av(:,:,:), aw(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: volp(:, :, :)

        REAL(realk) :: areae, areaw, arean, areas, areat, areab
        REAL(realk) :: vole, volw, voln, vols, volt, volb
        REAL(realk) :: mask

        ! Calculate cell volumes
        CALL get_field(au_f, "AU")
        CALL get_field(av_f, "AV")
        CALL get_field(aw_f, "AW")
        CALL set_field("VOLP")
        CALL this%calc_volp(this%bzelltyp, au_f, av_f, aw_f, this%icells, &
            this%icellspointer, this%xpsw, this%nvecs)


        ! 3-D fields used in the pressure solver
        CALL set_field("GSAW")
        CALL set_field("GSAE")
        CALL set_field("GSAS")
        CALL set_field("GSAN")
        CALL set_field("GSAB")
        CALL set_field("GSAT")
        CALL set_field("GSAP")

        DO igr = 1, nmygrids
            igrid = mygrids(igr)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL get_fieldptr(au, "AU", igrid)
            CALL get_fieldptr(av, "AV", igrid)
            CALL get_fieldptr(aw, "AW", igrid)

            CALL get_fieldptr(volp, "VOLP", igrid)

            CALL get_fieldptr(gsaw, "GSAW", igrid)
            CALL get_fieldptr(gsae, "GSAE", igrid)
            CALL get_fieldptr(gsas, "GSAS", igrid)
            CALL get_fieldptr(gsan, "GSAN", igrid)
            CALL get_fieldptr(gsab, "GSAB", igrid)
            CALL get_fieldptr(gsat, "GSAT", igrid)

            CALL get_fieldptr(gsap, "GSAP", igrid)

            ! Compute off-diagonal coefficients
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ! Open area of the cell faces
                        areae = au(k, j, i) * ddy(j) * ddz(k)
                        areaw = au(k, j, i-1) * ddy(j) * ddz(k)
                        arean = av(k, j, i) * ddx(i) * ddz(k)
                        areas = av(k, j-1, i) * ddx(i) * ddz(k)
                        areat = aw(k, j, i) * ddx(i) * ddy(j)
                        areab = aw(k-1, j, i) * ddx(i) * ddy(j)

                        ! Compute volume of the momentum cells
                        vole = (1/2)*(volp(k, j, i)+volp(k, j, i+1) &
                            +(1-CEILING(au(k, j, i-1))*volp(k, j, i) &
                            +(1-CEILING(au(k, j, i+1))*volp(k, j, i+1))

                        volw = (1/2)*(volp(k, j, i-1)+volp(k, j, i) &
                            +(1-CEILING(au(k, j, i-2))*volp(k, j, i-1) &
                            +(1-CEILING(au(k, j, i))*volp(k, j, i))

                        voln = (1/2)*(volp(k, j, i)+volp(k, j+1, i) &
                            +(1-CEILING(av(k, j-1, i))*volp(k, j, i) &
                            +(1-CEILING(av(k, j+1, i))*volp(k, j+1, i))

                        vols = (1/2)*(volp(k, j-1, i)+volp(k, j, i) &
                            +(1-CEILING(av(k, j-2, i))*volp(k, j-1, i) &
                            +(1-CEILING(av(k, j, i))*volp(k, j, i))

                        volt = (1/2)*(volp(k, j, i)+volp(k+1, j, i) &
                            +(1-CEILING(aw(k-1, j, i))*volp(k, j, i) &
                            +(1-CEILING(aw(k+1, j, i))*volp(k+1, j, i))

                        vols = (1/2)*(volp(k-1, j, i)+volp(k, j, i) &
                            +(1-CEILING(aw(k-2, j, i))*volp(k-1, j, i) &
                            +(1-CEILING(aw(k, j, i))*volp(k, j, i))

                        ! Only consider cells for which VOLP > 0
                        mask = CEILING(MIN(1.0, volp(k, j, i)/(ddx(i)*ddy(j)*ddz(k))))

                        gsae(k, j, i) = mask*divide0(areae, vole)
                        gsaw(k, j, i) = mask*divide0(areaw, volw)
                        gsan(k, j, i) = mask*divide0(arean, voln)
                        gsas(k, j, i) = mask*divide0(areas, vols)
                        gsat(k, j, i) = mask*divide0(areat, volt)
                        gsab(k, j, i) = mask*divide0(areab, volb)
                    END DO
                END DO
            END DO

            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            ! Front/West
            IF (nfro == 2 .OR. nfro == 5 .OR. nfro == 6 .OR. nfro == 19) THEN
                gsaw(k, j, 3) = 0.0
            END IF

            ! Back/East
            IF (nbac == 2 .OR. nbac == 5 .OR. nbac == 6) THEN
                gsae(k, j, ii-2) = 0.0
            END IF

            ! Right/South
            IF (nrgt == 2 .OR. nrgt == 5 .OR. nrgt == 6 .OR. nrgt == 19) THEN
                gsas(k, 3, i) = 0.0
            END IF

            ! Left/North
            IF (nlft == 2 .OR. nlft == 5 .OR. nlft == 6) THEN
                gsan(k, jj-2, i) = 0.0
            END IF

            ! Bottom
            IF (nbot == 2 .OR. nbot == 5 .OR. nbot == 6 .OR. nbot == 19) THEN
                gsab(3, j, i) = 0.0
            END IF

            ! Top
            IF (ntop == 2 .OR. ntop == 5 .OR. ntop == 6) THEN
                gsat(kk-2, j, i) = 0.0
            END IF
           
            ! Compute diagonal coefficient
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        gsap(k, j, i) = -gsae(k, j, i)-gsaw(k, j, i)-gsan(k, j, i) &
                            -gsas(k, j, i)-gsat(k, j, i)-gsab(k, j, i)
                    END DO
                END DO
            END DO
        END DO
    END SUBROUTINE giteig

    SUBROUTINE calc_nvecs(this, bzelltyp, au, av, aw, icells, icellspointer, &
            xpsw, nvecs, arealist)
        ! Subroutine arguments
        CLASS(gc_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(in) :: bzelltyp(*)
        TYPE(field_t), INTENT(in) :: au
        TYPE(field_t), INTENT(in) :: av
        TYPE(field_t), INTENT(in) :: aw
        INTEGER(intk), INTENT(in) :: icells(:)
        INTEGER(intk), INTENT(in) :: icellspointer(:)
        REAL(realk), INTENT(in) :: xpsw(:, :)
        REAL(realk), INTENT(out) :: nvecs(:, :)
        REAL(realk), INTENT(out) :: arealist(:)

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ipp, ncells
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            ipp = icellspointer(igrid)
            ncells = icells(igrid)

            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL this%calc_nvecs_grid(kk, jj, ii, ddx, ddy, ddz, &
                bzelltyp(ip3), au%arr(ip3), av%arr(ip3), aw%arr(ip3), &
                ncells, xpsw(:, ipp:ipp+ncells-1), &
                nvecs(:, ipp:ipp+ncells-1), arealist(ipp:ipp+ncells-1))
        END DO
    END SUBROUTINE calc_nvecs


    SUBROUTINE calc_nvecs_grid(kk, jj, ii, ddx, ddy, ddz, bzelltyp, &
            au, av, aw, ncells, xpsw, nvecs, arealist)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        REAL(realk), INTENT(in) :: au(kk, jj, ii)
        REAL(realk), INTENT(in) :: av(kk, jj, ii)
        REAL(realk), INTENT(in) :: aw(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: ncells
        REAL(realk), INTENT(in) :: xpsw(:, :)
        REAL(realk), INTENT(out) :: nvecs(:, :)
        REAL(realk), INTENT(out) :: arealist(:)

        ! Local variables
        INTEGER(intk) :: k, j, i, icell
        REAL(realk) :: nx, ny, nz, area, cartarea

        ! Corresponding to indices in calcauavaw_mod.F90 where this
        ! code initially came from
        DO i = 2, ii
            DO j = 2, jj
                DO k = 2, kk
                    IF (bzelltyp(k, j, i) >= 0) CYCLE
                    icell = -bzelltyp(k, j, i)

                    area = SQRT( &
                        (ddz(k)*ddy(j)*(au(k, j, i)-au(k, j, i-1)))**2 &
                        + (ddx(i)*ddz(k)*(av(k, j, i)-av(k, j-1, i)))**2 &
                        + (ddx(i)*ddy(j)*(aw(k, j, i)-aw(k-1, j, i)))**2)

                    ! From calcauavaw
                    ! Referenzflaeche: Mittel der kartesischen
                    ! Flaechen aus allen drei Richtungen.
                    cartarea = (1.0/3.0)*(ddx(i)*ddy(j) + &
                        ddx(i)*ddz(k) + ddz(k)*ddy(j))

                    ! Kann vorkommen, dass Flaeche zu klein,
                    ! dann muss Notfallnormale berechnet werden.
                    IF (area/cartarea < maccur**2) THEN
                        nx = 0.0
                        ny = 0.0
                        nz = 0.0
                    ELSE
                        nx = (au(k, j, i)-au(k, j, i-1))*ddz(k)*ddy(j)/area
                        ny = (av(k, j, i)-av(k, j-1, i))*ddx(i)*ddz(k)/area
                        nz = (aw(k, j, i)-aw(k-1, j, i))*ddx(i)*ddy(j)/area
                    END IF

                    arealist(icell) = area
                    nvecs(1, icell) = nx
                    nvecs(2, icell) = ny
                    nvecs(3, icell) = nz
                    nvecs(4, icell) = nx*xpsw(1, icell) &
                        + ny*xpsw(2, icell) &
                        + nz*xpsw(3, icell)
                END DO
            END DO
        END DO
    END SUBROUTINE calc_nvecs_grid


    SUBROUTINE calc_volp(this, bzelltyp, au, av, aw, icells, icellspointer, &
            xpsw, nvecs)
        ! Subroutine arguments
        CLASS(gc_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(in) :: bzelltyp(*)
        TYPE(field_t), INTENT(in) :: au
        TYPE(field_t), INTENT(in) :: av
        TYPE(field_t), INTENT(in) :: aw
        INTEGER(intk), INTENT(in) :: icells(:)
        INTEGER(intk), INTENT(in) :: icellspointer(:)
        REAL(realk), INTENT(in) :: xpsw(:, :)
        REAL(realk), INTENT(in) :: nvecs(:, :)

        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ipp, ncells
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)
        REAL(realk), POINTER, CONTIGUOUS :: volp(:,:,:)

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            ipp = icellspointer(igrid)
            ncells = icells(igrid)

            CALL get_fieldptr(x, "X", igrid)
            CALL get_fieldptr(y, "Y", igrid)
            CALL get_fieldptr(z, "Z", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)
            CALL get_fieldptr(volp, "VOLP", igrid)

            CALL this%calc_volp_grid(kk, jj, ii, x, y, z, ddx, ddy, ddz, &
                bzelltyp(ip3), au%arr(ip3), av%arr(ip3), aw%arr(ip3), &
                ncells, xpsw(:, ipp:ipp+ncells-1), nvecs(:, ipp:ipp+ncells-1) &
                volp)
        END DO
    END SUBROUTINE calc_volp

    SUBROUTINE calc_volp_grid(kk, jj, ii, x, y, z, ddx, ddy, ddz, bzelltyp, &
            au, av, aw, ncells, xpsw, nvecs, volp)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: bzelltyp(kk, jj, ii)
        REAL(realk), INTENT(in) :: au(kk, jj, ii)
        REAL(realk), INTENT(in) :: av(kk, jj, ii)
        REAL(realk), INTENT(in) :: aw(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: ncells
        REAL(realk), INTENT(in) :: xpsw(:, :)
        REAL(realk), INTENT(in) :: nvecs(:, :)
        REAL(realk), INTENT(out) :: volp(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i, icell
        REAL(realk) :: nx, ny, nz, dpl, vol

        ! Corresponding to indices in calcauavaw_mod.F90 where this
        ! code initially came from
        DO i = 2, ii
            DO j = 2, jj
                DO k = 2, kk
                    IF (bzelltyp(k, j, i) == 1) THEN
                        volp(k, j, i) = ddx(i)*ddy(j)*ddz(k)

                    ELSE IF (bzelltyp(k, j, i) == 0) THEN
                        volp(k, j, i) = 0

                    ELSE IF (bzelltyp(k, j, i) < 0) THEN
                        icell = -bzelltyp(k, j, i)

                        nx = nvecs(1, icell)
                        ny = nvecs(2, icell)
                        nz = nvecs(3, icell)
                        dpl = nvecs(4, icell)

                        CALL calc_volp_cell(x(i), y(j), z(k), &
                            ddx(i), ddy(j), ddz(k), &
                            nx, ny, nz, dpl, vol)

                        volp(k, j, i) = vol
                    ELSE
                        CYCLE
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE calc_volp_grid

    SUBROUTINE calc_volp_cell(x, y, z, ddx, ddy, ddz, nx, ny, nz, dpl, vol)
        ! Subroutine arguments
        REAL(realk), INTENT(in) :: x, y, z
        REAL(realk), INTENT(in) :: ddx, ddy, ddz
        REAL(realk), INTENT(in) :: nx, ny, nz
        REAL(realk), INTENT(in) :: dpl
        REAL(realk), INTENT(out) :: vol

        ! Local variables
        REAL(realk) :: d0
        REAL(realk) :: n1, n2, n3, d
        REAL(realk) :: s1, s2, s3
        REAL(realk) :: fac
        REAL(realk) :: V

        ! Reference for cube-plane intersection volume
        !   Lehmann, M. & Gekle, S.: Analytic Solution to the Piecewise Linear
        !   Interface Construction Problem and its Application in Curvature
        !   Calculation for Volume-of-Fluid Simulation Codes (2021)
        !   https://arxiv.org/pdf/2006.12838.pdf

        ! Move origin to cell centre
        d0 = dpl + nx*x + ny*y + nz*z

        ! Rescale the cell to side length 1
        nx = nx*ddx/SQRT(ddx**2+ddy**2+ddz**2)
        ny = ny*ddy/SQRT(ddx**2+ddy**2+ddz**2)
        nz = nz*ddz/SQRT(ddx**2+ddy**2+ddz**2)
        d0 = d0/SQRT(ddx**2+ddy**2+ddz**2)

        IF (d0<-(ABS(nx)+ABS(ny)+ABS(nz))/2) THEN
            vol = 0

        ELSE IF (d0>(ABS(nx)+ABS(ny)+ABS(nz))/2) THEN
            vol = ddx*ddy*ddz

        ELSE

            ! Orient coordinate system such that 0 <= n1 <= n2 <= n3 <= 1
            n1 = MIN(ABS(nx), ABS(ny), ABS(nz))
            n3 = MAX(ABS(nx), ABS(ny), ABS(nz))
            n2 = ABS(nx)+ABS(ny)+ABS(nz)-n1-n3

            ! Move origin to corner
            d = (n1+n2+n3)/2-ABS(d0)

            ! Compute volume
            fac= 1/(6*n1*n2*n3)

            IF (d<n1) THEN
                V = fac*d**3

            ELSE IF (n1<=d .AND. d<=n2) THEN
                V = fac*(d**3-(d-n1)**3)

            ELSE IF (n2<=d .AND. d<=MIN(n1+n2,n3)) THEN
                V = fac*(d**3-(d-n1)**3-(d-n2)**3)

            ELSE IF (n3<=d) THEN
                V = fac*(d**3-(d-n1)**3-(d-n2)**3-(d-n3)**3)

            ELSE IF (MIN(n1+n2,n3)<=d .AND. d<=n3) THEN
                V = fac*6*n1*n2*(d-(n1+n2)/2)

            ELSE
                V = 0

            END IF

            vol = ddx*ddy*ddz*(SIGN(d0)*(1/2-V)+1/2)

        END IF

    END SUBROUTINE

    SUBROUTINE divcal(this, div, u, v, w, fak, ctyp)
        ! Subroutine arguments
        CLASS(gc_t), INTENT(inout) :: this
        TYPE(field_t), INTENT(inout) :: div
        TYPE(field_t), INTENT(in) :: u
        TYPE(field_t), INTENT(in) :: v
        TYPE(field_t), INTENT(in) :: w
        REAL(realk), INTENT(in) :: fak
        CHARACTER(len=1), INTENT(in), OPTIONAL :: ctyp

        ! Local variables
        INTEGER(intk) :: i, igrid, ip3
        INTEGER(intk) :: kk, jj, ii
        LOGICAL :: use_sdiv
        TYPE(field_t), POINTER :: rddx_f, rddy_f, rddz_f
        TYPE(field_t), POINTER :: sdiv_f, bp_f
        REAL(realk), CONTIGUOUS, POINTER :: rddx(:), rddy(:), rddz(:)
        REAL(realk), CONTIGUOUS, POINTER :: sdiv(:, :, :), bp(:, :, :)

        ! For safety
        NULLIFY(bp_f)
        NULLIFY(sdiv_f)
        NULLIFY(sdiv)

        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")
        CALL get_field(bp_f, "BP")

        use_sdiv = .TRUE.
        IF (PRESENT(ctyp)) THEN
            IF (ctyp == "W") THEN
                use_sdiv = .FALSE.
            END IF
        END IF

        IF (use_sdiv) THEN
            CALL get_field(sdiv_f, "SDIV")
        END IF

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL rddx_f%get_ptr(rddx, igrid)
            CALL rddy_f%get_ptr(rddy, igrid)
            CALL rddz_f%get_ptr(rddz, igrid)
            CALL bp_f%get_ptr(bp, igrid)

            IF (use_sdiv) CALL sdiv_f%get_ptr(sdiv, igrid)

            CALL this%divcal_grid(kk, jj, ii, fak, div%arr(ip3), u%arr(ip3), &
                v%arr(ip3), w%arr(ip3), rddx, rddy, rddz, bp, sdiv)
        END DO
    END SUBROUTINE divcal

END MODULE gc_mod
