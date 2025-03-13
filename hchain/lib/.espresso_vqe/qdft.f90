#define ZERO (0.D0, 0.D0)

SUBROUTINE qdft(h_psi_ptr, npw, npwx, nvec, nvecx, evc)
    USE util_param,       ONLY : DP, stdout
    USE noncollin_module, ONLY : npol
    USE mp_bands_util,    ONLY : gstart

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx
    COMPLEX(DP), INTENT(INOUT) :: evc(npwx,nvec)
    INTEGER :: i, j, k, ierr
    COMPLEX(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:)
    EXTERNAL h_psi_ptr

    INTERFACE
        SUBROUTINE vqe(npw, npwx, nvec, nvecx, npol, psi, hpsi) BIND(C)
            USE ISO_C_BINDING
            INTEGER(C_INT), INTENT(IN) :: npw, npwx, nvec, nvecx, npol
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: & 
                psi(npw*npol,nvec), hpsi(npw*npol,nvec)
        END SUBROUTINE vqe
    END INTERFACE

    ALLOCATE(hpsi(npwx, nvecx), STAT=ierr)
    IF(ierr /= 0) CALL errore('qdft', 'cannot allocate hpsi', ABS(ierr))
    hpsi = ZERO

    ALLOCATE(psi(npwx, nvecx), STAT=ierr)
    IF(ierr /= 0) CALL errore('qdft', 'cannot allocate psi', ABS(ierr))
    psi = ZERO

    DO k=1,nvec
        psi(1,k) = evc(1,k)
        IF(gstart == 2) psi(1,k) = CMPLX(DBLE(psi(1,k)), 0.D0, kind=DP)
        DO i=2,npwx
            psi(i,k) = evc(i,k)
        END DO
    END DO

    WRITE(*, *) "qdft.f90", npw, nvec, npol
    CALL h_psi_ptr(npwx, npw, nvec, psi, hpsi)
    CALL vqe(npw, npwx, nvec, nvecx, npol, psi, hpsi)

    DEALLOCATE(psi)
    DEALLOCATE(hpsi)
END SUBROUTINE qdft
