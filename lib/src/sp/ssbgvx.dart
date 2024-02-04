      void ssbgvx(JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO ) {

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
      double               AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, TEST, UPPER, VALEIG, WANTZ;
      String             ORDER, VECT;
      int                I, IINFO, INDD, INDE, INDEE, INDISP, INDIWO, INDWRK, ITMP1, J, JJ, NSPLIT;
      double               TMP1;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMV, SLACPY, SPBSTF, SSBGST, SSBTRD, SSTEBZ, SSTEIN, SSTEQR, SSTERF, SSWAP, XERBLA
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
         xerbla('SSBGVX', -INFO );
         return;
      }

      // Quick return if possible

      M = 0;
      if (N == 0) return;

      // Form a split Cholesky factorization of B.

      spbstf(UPLO, N, KB, BB, LDBB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         return;
      }

      // Transform problem to standard eigenvalue problem.

      ssbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q, LDQ, WORK, IINFO );

      // Reduce symmetric band matrix to tridiagonal form.

      INDD = 1;
      INDE = INDD + N;
      INDWRK = INDE + N;
      if ( WANTZ ) {
         VECT = 'U';
      } else {
         VECT = 'N';
      }
      ssbtrd(VECT, UPLO, N, KA, AB, LDAB, WORK( INDD ), WORK( INDE ), Q, LDQ, WORK( INDWRK ), IINFO );

      // If all eigenvalues are desired and ABSTOL is less than or equal
      // to zero, then call SSTERF or SSTEQR.  If this fails for some
      // eigenvalue, then try SSTEBZ.

      TEST = false;
      if ( INDEIG ) {
         if ( IL == 1 && IU == N ) {
            TEST = true;
         }
      }
      if ( ( ALLEIG || TEST ) && ( ABSTOL <= ZERO ) ) {
         scopy(N, WORK( INDD ), 1, W, 1 );
         INDEE = INDWRK + 2*N;
         scopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
         if ( !WANTZ ) {
            ssterf(N, W, WORK( INDEE ), INFO );
         } else {
            slacpy('A', N, N, Q, LDQ, Z, LDZ );
            ssteqr(JOBZ, N, W, WORK( INDEE ), Z, LDZ, WORK( INDWRK ), INFO );
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
      // call SSTEIN.

      if ( WANTZ ) {
         ORDER = 'B';
      } else {
         ORDER = 'E';
      }
      INDISP = 1 + N;
      INDIWO = INDISP + N;
      sstebz(RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, WORK( INDD ), WORK( INDE ), M, NSPLIT, W, IWORK( 1 ), IWORK( INDISP ), WORK( INDWRK ), IWORK( INDIWO ), INFO );

      if ( WANTZ ) {
         sstein(N, WORK( INDD ), WORK( INDE ), M, W, IWORK( 1 ), IWORK( INDISP ), Z, LDZ, WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO );

         // Apply transformation matrix used in reduction to tridiagonal
         // form to eigenvectors returned by SSTEIN.

         for (J = 1; J <= M; J++) { // 20
            scopy(N, Z( 1, J ), 1, WORK( 1 ), 1 );
            sgemv('N', N, N, ONE, Q, LDQ, WORK, 1, ZERO, Z( 1, J ), 1 );
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
               sswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
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