      SUBROUTINE CHETRD_HB2ST( STAGE1, VECT, UPLO, N, KD, AB, LDAB, D, E, HOUS, LHOUS, WORK, LWORK, INFO )


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
      COMPLEX            AB( LDAB, * ), HOUS( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               RZERO
      COMPLEX            ZERO, ONE
      const              RZERO = 0.0E+0, ZERO = ( 0.0E+0, 0.0E+0 ), ONE  = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WANTQ, UPPER, AFTERS1;
      int                I, M, K, IB, SWEEPID, MYID, SHIFT, STT, ST, ED, STIND, EDIND, BLKLASTIND, COLPT, THED, STEPERCOL, GRSIZ, THGRSIZ, THGRNB, THGRID, NBTILES, TTYPE, TID, NTHREADS, ABDPOS, ABOFDPOS, DPOS, OFDPOS, AWPOS, INDA, INDW, APOS, SIZEA, LDA, INDV, INDTAU, SICEV, SIZETAU, LDV, LHMIN, LWMIN;
      REAL               ABSTMP
      COMPLEX            TMP
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHB2ST_KERNELS, CLACPY, CLASET, XERBLA
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

      IB       = ILAENV2STAGE( 2, 'CHETRD_HB2ST', VECT, N, KD, -1, -1 )
      if ( N.EQ.0 .OR. KD.LE.1 ) {
         LHMIN = 1
         LWMIN = 1
      } else {
         LHMIN = ILAENV2STAGE( 3, 'CHETRD_HB2ST', VECT, N, KD, IB, -1 )
         LWMIN = ILAENV2STAGE( 4, 'CHETRD_HB2ST', VECT, N, KD, IB, -1 )
      }

      if ( .NOT.AFTERS1 .AND. .NOT.LSAME( STAGE1, 'N' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( VECT, 'N' ) ) {
         INFO = -2
      } else if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( KD.LT.0 ) {
         INFO = -5
      } else if ( LDAB.LT.(KD+1) ) {
         INFO = -7
      } else if ( LHOUS.LT.LHMIN .AND. .NOT.LQUERY ) {
         INFO = -11
      } else if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
         INFO = -13
      }

      if ( INFO.EQ.0 ) {
         HOUS( 1 ) = SROUNDUP_LWORK( LHMIN )
         WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      }

      if ( INFO.NE.0 ) {
         xerbla('CHETRD_HB2ST', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
          HOUS( 1 ) = 1
          WORK( 1 ) = 1
          RETURN
      }

      // Determine pointer position

      LDV      = KD + IB
      SIZETAU  = 2 * N
      SICEV    = 2 * N
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
      // complex because D is double and the imaginary part should be 0)
      // and store it in D. A sequential code here is better or
      // in a parallel environment it might need two cores for D and E

      if ( KD.EQ.0 ) {
          for (I = 1; I <= N; I++) { // 30
              D( I ) = REAL( AB( ABDPOS, I ) )
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

      if ( KD.EQ.1 ) {
          for (I = 1; I <= N; I++) { // 50
              D( I ) = REAL( AB( ABDPOS, I ) )
          } // 50

          // make off-diagonal elements real and copy them to E

          if ( UPPER ) {
              for (I = 1; I <= N - 1; I++) { // 60
                  TMP = AB( ABOFDPOS, I+1 )
                  ABSTMP = ABS( TMP )
                  AB( ABOFDPOS, I+1 ) = ABSTMP
                  E( I ) = ABSTMP
                  if ( ABSTMP.NE.RZERO ) {
                     TMP = TMP / ABSTMP
                  } else {
                     TMP = ONE
                  }
                  IF( I.LT.N-1 ) AB( ABOFDPOS, I+2 ) = AB( ABOFDPOS, I+2 )*TMP
                   // IF( WANTZ ) THEN
                      // CALL CSCAL( N, CONJG( TMP ), Q( 1, I+1 ), 1 )
                   // END IF
              } // 60
          } else {
              for (I = 1; I <= N - 1; I++) { // 70
                 TMP = AB( ABOFDPOS, I )
                 ABSTMP = ABS( TMP )
                 AB( ABOFDPOS, I ) = ABSTMP
                 E( I ) = ABSTMP
                 if ( ABSTMP.NE.RZERO ) {
                    TMP = TMP / ABSTMP
                 } else {
                    TMP = ONE
                 }
                 IF( I.LT.N-1 ) AB( ABOFDPOS, I+1 ) = AB( ABOFDPOS, I+1 )*TMP
                  // IF( WANTQ ) THEN
                     // CALL CSCAL( N, TMP, Q( 1, I+1 ), 1 )
                  // END IF
              } // 70
          }

          HOUS( 1 ) = 1
          WORK( 1 ) = 1
          RETURN
      }

      // Main code start here.
      // Reduce the hermitian band of A to a tridiagonal matrix.

      THGRSIZ   = N
      GRSIZ     = 1
      SHIFT     = 3
      NBTILES   = CEILING( REAL(N)/REAL(KD) )
      STEPERCOL = CEILING( REAL(SHIFT)/REAL(GRSIZ) )
      THGRNB    = CEILING( REAL(N-1)/REAL(THGRSIZ) )

      clacpy("A", KD+1, N, AB, LDAB, WORK( APOS ), LDA );
      claset("A", KD,   N, ZERO, ZERO, WORK( AWPOS ), LDA );


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
              IF( STT.GT.ED ) EXIT
              for (M = 1; M <= STEPERCOL; M++) { // 120
                  ST = STT
                  for (SWEEPID = ST; SWEEPID <= ED; SWEEPID++) { // 130
                      for (K = 1; K <= GRSIZ; K++) { // 140
                          MYID  = (I-SWEEPID)*(STEPERCOL*GRSIZ) + (M-1)*GRSIZ + K
                          if ( MYID.EQ.1 ) {
                              TTYPE = 1
                          } else {
                              TTYPE = MOD( MYID, 2 ) + 2
                          }

                          if ( TTYPE.EQ.2 ) {
                              COLPT      = (MYID/2)*KD + SWEEPID
                              STIND      = COLPT-KD+1
                              EDIND      = MIN(COLPT,N)
                              BLKLASTIND = COLPT
                          } else {
                              COLPT      = ((MYID+1)/2)*KD + SWEEPID
                              STIND      = COLPT-KD+1
                              EDIND      = MIN(COLPT,N)
                              if ( ( STIND.GE.EDIND-1 ).AND. ( EDIND.EQ.N ) ) {
                                  BLKLASTIND = N
                              } else {
                                  BLKLASTIND = 0
                              }
                          }

                          // Call the kernel

#if defined(_OPENMP) && _OPENMP >= 201307
                          if ( TTYPE.NE.1 ) {
!$OMP TASK DEPEND(in:WORK(MYID+SHIFT-1))
!$OMP$     DEPEND(in:WORK(MYID-1))
!$OMP$     DEPEND(out:WORK(MYID))
                              TID      = OMP_GET_THREAD_NUM()
                              chb2st_kernels(UPLO, WANTQ, TTYPE, STIND, EDIND, SWEEPID, N, KD, IB, WORK ( INDA ), LDA, HOUS( INDV ), HOUS( INDTAU ), LDV, WORK( INDW + TID*KD ) );
!$OMP END TASK
                          } else {
!$OMP TASK DEPEND(in:WORK(MYID+SHIFT-1))
!$OMP$     DEPEND(out:WORK(MYID))
                              TID      = OMP_GET_THREAD_NUM()
                              chb2st_kernels(UPLO, WANTQ, TTYPE, STIND, EDIND, SWEEPID, N, KD, IB, WORK ( INDA ), LDA, HOUS( INDV ), HOUS( INDTAU ), LDV, WORK( INDW + TID*KD ) );
!$OMP END TASK
                          }
#else
                          chb2st_kernels(UPLO, WANTQ, TTYPE, STIND, EDIND, SWEEPID, N, KD, IB, WORK ( INDA ), LDA, HOUS( INDV ), HOUS( INDTAU ), LDV, WORK( INDW ) );
#endif
                          if ( BLKLASTIND.GE.(N-1) ) {
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
          D( I ) = REAL( WORK( DPOS+(I-1)*LDA ) )
      } // 150

      // Copy the off diagonal from A to E. Note that E is REAL thus only
      // the Real part is needed, the imaginary part should be zero.

      if ( UPPER ) {
          for (I = 1; I <= N-1; I++) { // 160
             E( I ) = REAL( WORK( OFDPOS+I*LDA ) )
          } // 160
      } else {
          for (I = 1; I <= N-1; I++) { // 170
             E( I ) = REAL( WORK( OFDPOS+(I-1)*LDA ) )
          } // 170
      }

      WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      RETURN

      // End of CHETRD_HB2ST

      }
