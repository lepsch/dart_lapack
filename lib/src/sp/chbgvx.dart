      void chbgvx(JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, IFAIL, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, KA, KB, LDAB, LDBB, LDQ, LDZ, M, N;
      double               ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      double               RWORK( * ), W( * );
      Complex            AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, TEST, UPPER, VALEIG, WANTZ;
      String             ORDER, VECT;
      int                I, IINFO, INDD, INDE, INDEE, INDISP, INDIWK, INDRWK, INDWRK, ITMP1, J, JJ, NSPLIT;
      double               TMP1;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMV, CHBGST, CHBTRD, CLACPY, CPBSTF, CSTEIN, CSTEQR, CSWAP, SCOPY, SSTEBZ, SSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      UPPER = lsame( UPLO, 'U' );
      ALLEIG = lsame( RANGE, 'A' );
      VALEIG = lsame( RANGE, 'V' );
      INDEIG = lsame( RANGE, 'I' );

      INFO = 0;
      if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( ALLEIG || VALEIG || INDEIG ) ) {
         INFO = -2;
      } else if ( !( UPPER || lsame( UPLO, 'L' ) ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( KA < 0 ) {
         INFO = -5;
      } else if ( KB < 0 || KB > KA ) {
         INFO = -6;
      } else if ( LDAB < KA+1 ) {
         INFO = -8;
      } else if ( LDBB < KB+1 ) {
         INFO = -10;
      } else if ( LDQ < 1 || ( WANTZ && LDQ < N ) ) {
         INFO = -12;
      } else {
         if ( VALEIG ) {
            if (N > 0 && VU <= VL) INFO = -14;
         } else if ( INDEIG ) {
            if ( IL < 1 || IL > max( 1, N ) ) {
               INFO = -15;
            } else if ( IU < min( N, IL ) || IU > N ) {
               INFO = -16;
            }
         }
      }
      if ( INFO == 0) {
         if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
            INFO = -21;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CHBGVX', -INFO );
         return;
      }

      // Quick return if possible

      M = 0;
      if (N == 0) return;

      // Form a split Cholesky factorization of B.

      cpbstf(UPLO, N, KB, BB, LDBB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         return;
      }

      // Transform problem to standard eigenvalue problem.

      chbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q, LDQ, WORK, RWORK, IINFO );

      // Solve the standard eigenvalue problem.
      // Reduce Hermitian band matrix to tridiagonal form.

      INDD = 1;
      INDE = INDD + N;
      INDRWK = INDE + N;
      INDWRK = 1;
      if ( WANTZ ) {
         VECT = 'U';
      } else {
         VECT = 'N';
      }
      chbtrd(VECT, UPLO, N, KA, AB, LDAB, RWORK( INDD ), RWORK( INDE ), Q, LDQ, WORK( INDWRK ), IINFO );

      // If all eigenvalues are desired and ABSTOL is less than or equal
      // to zero, then call SSTERF or CSTEQR.  If this fails for some
      // eigenvalue, then try SSTEBZ.

      TEST = false;
      if ( INDEIG ) {
         if ( IL == 1 && IU == N ) {
            TEST = true;
         }
      }
      if ( ( ALLEIG || TEST ) && ( ABSTOL <= ZERO ) ) {
         scopy(N, RWORK( INDD ), 1, W, 1 );
         INDEE = INDRWK + 2*N;
         scopy(N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 );
         if ( !WANTZ ) {
            ssterf(N, W, RWORK( INDEE ), INFO );
         } else {
            clacpy('A', N, N, Q, LDQ, Z, LDZ );
            csteqr(JOBZ, N, W, RWORK( INDEE ), Z, LDZ, RWORK( INDRWK ), INFO );
            if ( INFO == 0 ) {
               for (I = 1; I <= N; I++) { // 10
                  IFAIL[I] = 0;
               } // 10
            }
         }
         if ( INFO == 0 ) {
            M = N;
            GO TO 30;
         }
         INFO = 0;
      }

      // Otherwise, call SSTEBZ and, if eigenvectors are desired,
      // call CSTEIN.

      if ( WANTZ ) {
         ORDER = 'B';
      } else {
         ORDER = 'E';
      }
      INDISP = 1 + N;
      INDIWK = INDISP + N;
      sstebz(RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, RWORK( INDD ), RWORK( INDE ), M, NSPLIT, W, IWORK( 1 ), IWORK( INDISP ), RWORK( INDRWK ), IWORK( INDIWK ), INFO );

      if ( WANTZ ) {
         cstein(N, RWORK( INDD ), RWORK( INDE ), M, W, IWORK( 1 ), IWORK( INDISP ), Z, LDZ, RWORK( INDRWK ), IWORK( INDIWK ), IFAIL, INFO );

         // Apply unitary matrix used in reduction to tridiagonal
         // form to eigenvectors returned by CSTEIN.

         for (J = 1; J <= M; J++) { // 20
            ccopy(N, Z( 1, J ), 1, WORK( 1 ), 1 );
            cgemv('N', N, N, CONE, Q, LDQ, WORK, 1, CZERO, Z( 1, J ), 1 );
         } // 20
      }

      } // 30

      // If eigenvalues are not in order, then sort them, along with
      // eigenvectors.

      if ( WANTZ ) {
         for (J = 1; J <= M - 1; J++) { // 50
            I = 0;
            TMP1 = W( J );
            for (JJ = J + 1; JJ <= M; JJ++) { // 40
               if ( W( JJ ) < TMP1 ) {
                  I = JJ;
                  TMP1 = W( JJ );
               }
            } // 40

            if ( I != 0 ) {
               ITMP1 = IWORK( 1 + I-1 );
               W[I] = W( J );
               IWORK[1 + I-1] = IWORK( 1 + J-1 );
               W[J] = TMP1;
               IWORK[1 + J-1] = ITMP1;
               cswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
               if ( INFO != 0 ) {
                  ITMP1 = IFAIL( I );
                  IFAIL[I] = IFAIL( J );
                  IFAIL[J] = ITMP1;
               }
            }
         } // 50
      }

      return;
      }