      SUBROUTINE SSYTRD_SB2ST( STAGE1, VECT, UPLO, N, KD, AB, LDAB, D, E, HOUS, LHOUS, WORK, LWORK, INFO )

#if defined(_OPENMP)
      use omp_lib
#endif

      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             STAGE1, UPLO, VECT;
      int                N, KD, LDAB, LHOUS, LWORK, INFO;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * )
      REAL               AB( LDAB, * ), HOUS( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               RZERO
      REAL               ZERO, ONE
      const              RZERO = 0.0E+0, ZERO = 0.0E+0, ONE  = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WANTQ, UPPER, AFTERS1;
      int                I, M, K, IB, SWEEPID, MYID, SHIFT, STT, ST, ED, STIND, EDIND, BLKLASTIND, COLPT, THED, STEPERCOL, GRSIZ, THGRSIZ, THGRNB, THGRID, NBTILES, TTYPE, TID, NTHREADS, ABDPOS, ABOFDPOS, DPOS, OFDPOS, AWPOS, INDA, INDW, APOS, SIZEA, LDA, INDV, INDTAU, SISEV, SIZETAU, LDV, LHMIN, LWMIN;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSB2ST_KERNELS, SLACPY, SLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, MAX, CEILING, REAL
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV2STAGE, SROUNDUP_LWORK
      // ..
      // .. Executable Statements ..

      // Determine the minimal workspace size required.
      // Test the input parameters

      INFO    = 0
      AFTERS1 = LSAME( STAGE1, 'Y' )
      WANTQ   = LSAME( VECT, 'V' )
      UPPER   = LSAME( UPLO, 'U' )
      LQUERY  = ( LWORK.EQ.-1 ) .OR. ( LHOUS.EQ.-1 )

      // Determine the block size, the workspace size and the hous size.

      IB       = ILAENV2STAGE( 2, 'SSYTRD_SB2ST', VECT, N, KD, -1, -1 )
      IF( N.EQ.0 .OR. KD.LE.1 ) THEN
         LHMIN = 1
         LWMIN = 1
      ELSE
         LHMIN = ILAENV2STAGE( 3, 'SSYTRD_SB2ST', VECT, N, KD, IB, -1 )
         LWMIN = ILAENV2STAGE( 4, 'SSYTRD_SB2ST', VECT, N, KD, IB, -1 )
      END IF

      IF( .NOT.AFTERS1 .AND. .NOT.LSAME( STAGE1, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( VECT, 'N' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( KD.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.(KD+1) ) THEN
         INFO = -7
      ELSE IF( LHOUS.LT.LHMIN .AND. .NOT.LQUERY ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -13
      END IF

      IF( INFO.EQ.0 ) THEN
         HOUS( 1 ) = SROUNDUP_LWORK( LHMIN )
         WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYTRD_SB2ST', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) THEN
          HOUS( 1 ) = 1
          WORK( 1 ) = 1
          RETURN
      END IF

      // Determine pointer position

      LDV      = KD + IB
      SIZETAU  = 2 * N
      SISEV    = 2 * N
      INDTAU   = 1
      INDV     = INDTAU + SIZETAU
      LDA      = 2 * KD + 1
      SIZEA    = LDA * N
      INDA     = 1
      INDW     = INDA + SIZEA
      NTHREADS = 1
      TID      = 0

      IF( UPPER ) THEN
          APOS     = INDA + KD
          AWPOS    = INDA
          DPOS     = APOS + KD
          OFDPOS   = DPOS - 1
          ABDPOS   = KD + 1
          ABOFDPOS = KD
      ELSE
          APOS     = INDA
          AWPOS    = INDA + KD + 1
          DPOS     = APOS
          OFDPOS   = DPOS + 1
          ABDPOS   = 1
          ABOFDPOS = 2

      ENDIF

      // Case KD=0:
      // The matrix is diagonal. We just copy it (convert to "real" for
      // real because D is double and the imaginary part should be 0)
      // and store it in D. A sequential code here is better or
      // in a parallel environment it might need two cores for D and E

      IF( KD.EQ.0 ) THEN
          DO 30 I = 1, N
              D( I ) = ( AB( ABDPOS, I ) )
   30     CONTINUE
          DO 40 I = 1, N-1
              E( I ) = RZERO
   40     CONTINUE

          HOUS( 1 ) = 1
          WORK( 1 ) = 1
          RETURN
      END IF

      // Case KD=1:
      // The matrix is already Tridiagonal. We have to make diagonal
      // and offdiagonal elements real, and store them in D and E.
      // For that, for real precision just copy the diag and offdiag
     t // o D and E while for the COMPLEX case the bulge chasing is
      // performed to convert the hermetian tridiagonal to symmetric
     t // ridiagonal. A simpler conversion formula might be used, but then
      // updating the Q matrix will be required and based if Q is generated
      // or not this might complicate the story.

      IF( KD.EQ.1 ) THEN
          DO 50 I = 1, N
              D( I ) = ( AB( ABDPOS, I ) )
   50     CONTINUE

          IF( UPPER ) THEN
              DO 60 I = 1, N-1
                 E( I ) = ( AB( ABOFDPOS, I+1 ) )
   60         CONTINUE
          ELSE
              DO 70 I = 1, N-1
                 E( I ) = ( AB( ABOFDPOS, I ) )
   70         CONTINUE
          ENDIF

          HOUS( 1 ) = 1
          WORK( 1 ) = 1
          RETURN
      END IF

      // Main code start here.
      // Reduce the symmetric band of A to a tridiagonal matrix.

      THGRSIZ   = N
      GRSIZ     = 1
      SHIFT     = 3
      NBTILES   = CEILING( REAL(N)/REAL(KD) )
      STEPERCOL = CEILING( REAL(SHIFT)/REAL(GRSIZ) )
      THGRNB    = CEILING( REAL(N-1)/REAL(THGRSIZ) )

      CALL SLACPY( "A", KD+1, N, AB, LDAB, WORK( APOS ), LDA )
      CALL SLASET( "A", KD,   N, ZERO, ZERO, WORK( AWPOS ), LDA )


      // openMP parallelisation start here

#if defined(_OPENMP)
!$OMP PARALLEL PRIVATE( TID, THGRID, BLKLASTIND )
!$OMP$         PRIVATE( THED, I, M, K, ST, ED, STT, SWEEPID )
!$OMP$         PRIVATE( MYID, TTYPE, COLPT, STIND, EDIND )
!$OMP$         SHARED ( UPLO, WANTQ, INDV, INDTAU, HOUS, WORK)
!$OMP$         SHARED ( N, KD, IB, NBTILES, LDA, LDV, INDA )
!$OMP$         SHARED ( STEPERCOL, THGRNB, THGRSIZ, GRSIZ, SHIFT )
!$OMP MASTER
#endif

      // main bulge chasing loop

      DO 100 THGRID = 1, THGRNB
          STT  = (THGRID-1)*THGRSIZ+1
          THED = MIN( (STT + THGRSIZ -1), (N-1))
          DO 110 I = STT, N-1
              ED = MIN( I, THED )
              IF( STT.GT.ED ) EXIT
              DO 120 M = 1, STEPERCOL
                  ST = STT
                  DO 130 SWEEPID = ST, ED
                      DO 140 K = 1, GRSIZ
                          MYID  = (I-SWEEPID)*(STEPERCOL*GRSIZ) + (M-1)*GRSIZ + K
                          IF ( MYID.EQ.1 ) THEN
                              TTYPE = 1
                          ELSE
                              TTYPE = MOD( MYID, 2 ) + 2
                          ENDIF

                          IF( TTYPE.EQ.2 ) THEN
                              COLPT      = (MYID/2)*KD + SWEEPID
                              STIND      = COLPT-KD+1
                              EDIND      = MIN(COLPT,N)
                              BLKLASTIND = COLPT
                          ELSE
                              COLPT      = ((MYID+1)/2)*KD + SWEEPID
                              STIND      = COLPT-KD+1
                              EDIND      = MIN(COLPT,N)
                              IF( ( STIND.GE.EDIND-1 ).AND. ( EDIND.EQ.N ) ) THEN
                                  BLKLASTIND = N
                              ELSE
                                  BLKLASTIND = 0
                              ENDIF
                          ENDIF

                          // Call the kernel

#if defined(_OPENMP) && _OPENMP >= 201307
                          IF( TTYPE.NE.1 ) THEN
!$OMP TASK DEPEND(in:WORK(MYID+SHIFT-1))
!$OMP$     DEPEND(in:WORK(MYID-1))
!$OMP$     DEPEND(out:WORK(MYID))
                              TID      = OMP_GET_THREAD_NUM()
                              CALL SSB2ST_KERNELS( UPLO, WANTQ, TTYPE, STIND, EDIND, SWEEPID, N, KD, IB, WORK ( INDA ), LDA, HOUS( INDV ), HOUS( INDTAU ), LDV, WORK( INDW + TID*KD ) )
!$OMP END TASK
                          ELSE
!$OMP TASK DEPEND(in:WORK(MYID+SHIFT-1))
!$OMP$     DEPEND(out:WORK(MYID))
                              TID      = OMP_GET_THREAD_NUM()
                              CALL SSB2ST_KERNELS( UPLO, WANTQ, TTYPE, STIND, EDIND, SWEEPID, N, KD, IB, WORK ( INDA ), LDA, HOUS( INDV ), HOUS( INDTAU ), LDV, WORK( INDW + TID*KD ) )
!$OMP END TASK
                          ENDIF
#else
                          CALL SSB2ST_KERNELS( UPLO, WANTQ, TTYPE, STIND, EDIND, SWEEPID, N, KD, IB, WORK ( INDA ), LDA, HOUS( INDV ), HOUS( INDTAU ), LDV, WORK( INDW ) )
#endif
                          IF ( BLKLASTIND.GE.(N-1) ) THEN
                              STT = STT + 1
                              EXIT
                          ENDIF
  140                 CONTINUE
  130             CONTINUE
  120         CONTINUE
  110     CONTINUE
  100 CONTINUE

#if defined(_OPENMP)
!$OMP END MASTER
!$OMP END PARALLEL
#endif

      // Copy the diagonal from A to D. Note that D is REAL thus only
     t // he Real part is needed, the imaginary part should be zero.

      DO 150 I = 1, N
          D( I ) = ( WORK( DPOS+(I-1)*LDA ) )
  150 CONTINUE

      // Copy the off diagonal from A to E. Note that E is REAL thus only
     t // he Real part is needed, the imaginary part should be zero.

      IF( UPPER ) THEN
          DO 160 I = 1, N-1
             E( I ) = ( WORK( OFDPOS+I*LDA ) )
  160     CONTINUE
      ELSE
          DO 170 I = 1, N-1
             E( I ) = ( WORK( OFDPOS+(I-1)*LDA ) )
  170     CONTINUE
      ENDIF

      WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      RETURN

      // End of SSYTRD_SB2ST

      }
