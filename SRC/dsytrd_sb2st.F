      SUBROUTINE DSYTRD_SB2ST( STAGE1, VECT, UPLO, N, KD, AB, LDAB, D, E, HOUS, LHOUS, WORK, LWORK, INFO )

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
      double             D( * ), E( * );
      double             AB( LDAB, * ), HOUS( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             RZERO;
      double             ZERO, ONE;
      const              RZERO = 0.0D+0, ZERO = 0.0D+0, ONE  = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WANTQ, UPPER, AFTERS1;
      int                I, M, K, IB, SWEEPID, MYID, SHIFT, STT, ST, ED, STIND, EDIND, BLKLASTIND, COLPT, THED, STEPERCOL, GRSIZ, THGRSIZ, THGRNB, THGRID, NBTILES, TTYPE, TID, NTHREADS, ABDPOS, ABOFDPOS, DPOS, OFDPOS, AWPOS, INDA, INDW, APOS, SIZEA, LDA, INDV, INDTAU, SIDEV, SIZETAU, LDV, LHMIN, LWMIN;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSB2ST_KERNELS, DLACPY, DLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, MAX, CEILING, REAL
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      // EXTERNAL LSAME, ILAENV2STAGE
      // ..
      // .. Executable Statements ..

      // Determine the minimal workspace size required.
      // Test the input parameters

      INFO    = 0
      AFTERS1 = LSAME( STAGE1, 'Y' )
      WANTQ   = LSAME( VECT, 'V' )
      UPPER   = LSAME( UPLO, 'U' )
      LQUERY  = ( LWORK == -1 ) || ( LHOUS == -1 )

      // Determine the block size, the workspace size and the hous size.

      IB       = ILAENV2STAGE( 2, 'DSYTRD_SB2ST', VECT, N, KD, -1, -1 )
      if ( N == 0 || KD.LE.1 ) {
         LHMIN = 1
         LWMIN = 1
      } else {
         LHMIN = ILAENV2STAGE( 3, 'DSYTRD_SB2ST', VECT, N, KD, IB, -1 )
         LWMIN = ILAENV2STAGE( 4, 'DSYTRD_SB2ST', VECT, N, KD, IB, -1 )
      }

      if ( .NOT.AFTERS1 && .NOT.LSAME( STAGE1, 'N' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( VECT, 'N' ) ) {
         INFO = -2
      } else if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -3
      } else if ( N < 0 ) {
         INFO = -4
      } else if ( KD < 0 ) {
         INFO = -5
      } else if ( LDAB < (KD+1) ) {
         INFO = -7
      } else if ( LHOUS < LHMIN && .NOT.LQUERY ) {
         INFO = -11
      } else if ( LWORK < LWMIN && .NOT.LQUERY ) {
         INFO = -13
      }

      if ( INFO == 0 ) {
         HOUS( 1 ) = LHMIN
         WORK( 1 ) = LWMIN
      }

      if ( INFO != 0 ) {
         xerbla('DSYTRD_SB2ST', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N == 0 ) {
          HOUS( 1 ) = 1
          WORK( 1 ) = 1
          RETURN
      }

      // Determine pointer position

      LDV      = KD + IB
      SIZETAU  = 2 * N
      SIDEV    = 2 * N
      INDTAU   = 1
      INDV     = INDTAU + SIZETAU
      LDA      = 2 * KD + 1
      SIZEA    = LDA * N
      INDA     = 1
      INDW     = INDA + SIZEA
      NTHREADS = 1
      TID      = 0

      if ( UPPER ) {
          APOS     = INDA + KD
          AWPOS    = INDA
          DPOS     = APOS + KD
          OFDPOS   = DPOS - 1
          ABDPOS   = KD + 1
          ABOFDPOS = KD
      } else {
          APOS     = INDA
          AWPOS    = INDA + KD + 1
          DPOS     = APOS
          OFDPOS   = DPOS + 1
          ABDPOS   = 1
          ABOFDPOS = 2

      }

      // Case KD=0:
      // The matrix is diagonal. We just copy it (convert to "real" for
      // real because D is double and the imaginary part should be 0)
      // and store it in D. A sequential code here is better or
      // in a parallel environment it might need two cores for D and E

      if ( KD == 0 ) {
          for (I = 1; I <= N; I++) { // 30
              D( I ) = ( AB( ABDPOS, I ) )
          } // 30
          for (I = 1; I <= N-1; I++) { // 40
              E( I ) = RZERO
          } // 40

          HOUS( 1 ) = 1
          WORK( 1 ) = 1
          RETURN
      }

      // Case KD=1:
      // The matrix is already Tridiagonal. We have to make diagonal
      // and offdiagonal elements real, and store them in D and E.
      // For that, for real precision just copy the diag and offdiag
      // to D and E while for the COMPLEX case the bulge chasing is
      // performed to convert the hermetian tridiagonal to symmetric
      // tridiagonal. A simpler conversion formula might be used, but then
      // updating the Q matrix will be required and based if Q is generated
      // or not this might complicate the story.

      if ( KD == 1 ) {
          for (I = 1; I <= N; I++) { // 50
              D( I ) = ( AB( ABDPOS, I ) )
          } // 50

          if ( UPPER ) {
              for (I = 1; I <= N-1; I++) { // 60
                 E( I ) = ( AB( ABOFDPOS, I+1 ) )
              } // 60
          } else {
              for (I = 1; I <= N-1; I++) { // 70
                 E( I ) = ( AB( ABOFDPOS, I ) )
              } // 70
          }

          HOUS( 1 ) = 1
          WORK( 1 ) = 1
          RETURN
      }

      // Main code start here.
      // Reduce the symmetric band of A to a tridiagonal matrix.

      THGRSIZ   = N
      GRSIZ     = 1
      SHIFT     = 3
      NBTILES   = CEILING( REAL(N)/REAL(KD) )
      STEPERCOL = CEILING( REAL(SHIFT)/REAL(GRSIZ) )
      THGRNB    = CEILING( REAL(N-1)/REAL(THGRSIZ) )

      dlacpy("A", KD+1, N, AB, LDAB, WORK( APOS ), LDA );
      dlaset("A", KD,   N, ZERO, ZERO, WORK( AWPOS ), LDA );


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

      for (THGRID = 1; THGRID <= THGRNB; THGRID++) { // 100
          STT  = (THGRID-1)*THGRSIZ+1
          THED = MIN( (STT + THGRSIZ -1), (N-1))
          for (I = STT; I <= N-1; I++) { // 110
              ED = MIN( I, THED )
              if (STT > ED) EXIT;
              for (M = 1; M <= STEPERCOL; M++) { // 120
                  ST = STT
                  for (SWEEPID = ST; SWEEPID <= ED; SWEEPID++) { // 130
                      for (K = 1; K <= GRSIZ; K++) { // 140
                          MYID  = (I-SWEEPID)*(STEPERCOL*GRSIZ) + (M-1)*GRSIZ + K
                          if ( MYID == 1 ) {
                              TTYPE = 1
                          } else {
                              TTYPE = MOD( MYID, 2 ) + 2
                          }

                          if ( TTYPE == 2 ) {
                              COLPT      = (MYID/2)*KD + SWEEPID
                              STIND      = COLPT-KD+1
                              EDIND      = MIN(COLPT,N)
                              BLKLASTIND = COLPT
                          } else {
                              COLPT      = ((MYID+1)/2)*KD + SWEEPID
                              STIND      = COLPT-KD+1
                              EDIND      = MIN(COLPT,N)
                              if ( ( STIND >= EDIND-1 ) && ( EDIND == N ) ) {
                                  BLKLASTIND = N
                              } else {
                                  BLKLASTIND = 0
                              }
                          }

                          // Call the kernel

#if defined(_OPENMP) &&  _OPENMP >= 201307
                          if ( TTYPE != 1 ) {
!$OMP TASK DEPEND(in:WORK(MYID+SHIFT-1))
!$OMP$     DEPEND(in:WORK(MYID-1))
!$OMP$     DEPEND(out:WORK(MYID))
                              TID      = OMP_GET_THREAD_NUM()
                              dsb2st_kernels(UPLO, WANTQ, TTYPE, STIND, EDIND, SWEEPID, N, KD, IB, WORK ( INDA ), LDA, HOUS( INDV ), HOUS( INDTAU ), LDV, WORK( INDW + TID*KD ) );
!$OMP END TASK
                          } else {
!$OMP TASK DEPEND(in:WORK(MYID+SHIFT-1))
!$OMP$     DEPEND(out:WORK(MYID))
                              TID      = OMP_GET_THREAD_NUM()
                              dsb2st_kernels(UPLO, WANTQ, TTYPE, STIND, EDIND, SWEEPID, N, KD, IB, WORK ( INDA ), LDA, HOUS( INDV ), HOUS( INDTAU ), LDV, WORK( INDW + TID*KD ) );
!$OMP END TASK
                          }
#else
                          dsb2st_kernels(UPLO, WANTQ, TTYPE, STIND, EDIND, SWEEPID, N, KD, IB, WORK ( INDA ), LDA, HOUS( INDV ), HOUS( INDTAU ), LDV, WORK( INDW ) );
#endif
                          if ( BLKLASTIND >= (N-1) ) {
                              STT = STT + 1
                              EXIT
                          }
                      } // 140
                  } // 130
              } // 120
          } // 110
      } // 100

#if defined(_OPENMP)
!$OMP END MASTER
!$OMP END PARALLEL
#endif

      // Copy the diagonal from A to D. Note that D is REAL thus only
      // the Real part is needed, the imaginary part should be zero.

      for (I = 1; I <= N; I++) { // 150
          D( I ) = ( WORK( DPOS+(I-1)*LDA ) )
      } // 150

      // Copy the off diagonal from A to E. Note that E is REAL thus only
      // the Real part is needed, the imaginary part should be zero.

      if ( UPPER ) {
          for (I = 1; I <= N-1; I++) { // 160
             E( I ) = ( WORK( OFDPOS+I*LDA ) )
          } // 160
      } else {
          for (I = 1; I <= N-1; I++) { // 170
             E( I ) = ( WORK( OFDPOS+(I-1)*LDA ) )
          } // 170
      }

      WORK( 1 ) = LWMIN
      RETURN

      // End of DSYTRD_SB2ST

      }
