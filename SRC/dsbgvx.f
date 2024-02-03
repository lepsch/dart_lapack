      SUBROUTINE DSBGVX( JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, KA, KB, LDAB, LDBB, LDQ, LDZ, M, N;
      double             ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      double             AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, TEST, UPPER, VALEIG, WANTZ;
      String             ORDER, VECT;
      int                I, IINFO, INDD, INDE, INDEE, INDISP, INDIWO, INDWRK, ITMP1, J, JJ, NSPLIT;
      double             TMP1;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMV, DLACPY, DPBSTF, DSBGST, DSBTRD, DSTEBZ, DSTEIN, DSTEQR, DSTERF, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )

      INFO = 0
      if ( .NOT.( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( ALLEIG || VALEIG || INDEIG ) ) {
         INFO = -2
      } else if ( .NOT.( UPPER || LSAME( UPLO, 'L' ) ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( KA.LT.0 ) {
         INFO = -5
      } else if ( KB.LT.0 || KB.GT.KA ) {
         INFO = -6
      } else if ( LDAB.LT.KA+1 ) {
         INFO = -8
      } else if ( LDBB.LT.KB+1 ) {
         INFO = -10
      } else if ( LDQ.LT.1 || ( WANTZ && LDQ.LT.N ) ) {
         INFO = -12
      } else {
         if ( VALEIG ) {
            if (N.GT.0 && VU.LE.VL) INFO = -14;
         } else if ( INDEIG ) {
            if ( IL.LT.1 || IL.GT.MAX( 1, N ) ) {
               INFO = -15
            } else if ( IU.LT.MIN( N, IL ) || IU.GT.N ) {
               INFO = -16
            }
         }
      }
      if ( INFO == 0) {
         if ( LDZ.LT.1 || ( WANTZ && LDZ.LT.N ) ) {
            INFO = -21
         }
      }

      if ( INFO != 0 ) {
         xerbla('DSBGVX', -INFO );
         RETURN
      }

      // Quick return if possible

      M = 0
      if (N == 0) RETURN;

      // Form a split Cholesky factorization of B.

      dpbstf(UPLO, N, KB, BB, LDBB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO
         RETURN
      }

      // Transform problem to standard eigenvalue problem.

      dsbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q, LDQ, WORK, IINFO );

      // Reduce symmetric band matrix to tridiagonal form.

      INDD = 1
      INDE = INDD + N
      INDWRK = INDE + N
      if ( WANTZ ) {
         VECT = 'U'
      } else {
         VECT = 'N'
      }
      dsbtrd(VECT, UPLO, N, KA, AB, LDAB, WORK( INDD ), WORK( INDE ), Q, LDQ, WORK( INDWRK ), IINFO );

      // If all eigenvalues are desired and ABSTOL is less than or equal
      // to zero, then call DSTERF or SSTEQR.  If this fails for some
      // eigenvalue, then try DSTEBZ.

      TEST = false;
      if ( INDEIG ) {
         if ( IL == 1 && IU == N ) {
            TEST = true;
         }
      }
      if ( ( ALLEIG || TEST ) && ( ABSTOL.LE.ZERO ) ) {
         dcopy(N, WORK( INDD ), 1, W, 1 );
         INDEE = INDWRK + 2*N
         dcopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
         if ( .NOT.WANTZ ) {
            dsterf(N, W, WORK( INDEE ), INFO );
         } else {
            dlacpy('A', N, N, Q, LDQ, Z, LDZ );
            dsteqr(JOBZ, N, W, WORK( INDEE ), Z, LDZ, WORK( INDWRK ), INFO );
            if ( INFO == 0 ) {
               for (I = 1; I <= N; I++) { // 10
                  IFAIL( I ) = 0
               } // 10
            }
         }
         if ( INFO == 0 ) {
            M = N
            GO TO 30
         }
         INFO = 0
      }

      // Otherwise, call DSTEBZ and, if eigenvectors are desired,
      // call DSTEIN.

      if ( WANTZ ) {
         ORDER = 'B'
      } else {
         ORDER = 'E'
      }
      INDISP = 1 + N
      INDIWO = INDISP + N
      dstebz(RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, WORK( INDD ), WORK( INDE ), M, NSPLIT, W, IWORK( 1 ), IWORK( INDISP ), WORK( INDWRK ), IWORK( INDIWO ), INFO );

      if ( WANTZ ) {
         dstein(N, WORK( INDD ), WORK( INDE ), M, W, IWORK( 1 ), IWORK( INDISP ), Z, LDZ, WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO );

         // Apply transformation matrix used in reduction to tridiagonal
         // form to eigenvectors returned by DSTEIN.

         for (J = 1; J <= M; J++) { // 20
            dcopy(N, Z( 1, J ), 1, WORK( 1 ), 1 );
            dgemv('N', N, N, ONE, Q, LDQ, WORK, 1, ZERO, Z( 1, J ), 1 );
         } // 20
      }

      } // 30

      // If eigenvalues are not in order, then sort them, along with
      // eigenvectors.

      if ( WANTZ ) {
         for (J = 1; J <= M - 1; J++) { // 50
            I = 0
            TMP1 = W( J )
            for (JJ = J + 1; JJ <= M; JJ++) { // 40
               if ( W( JJ ).LT.TMP1 ) {
                  I = JJ
                  TMP1 = W( JJ )
               }
            } // 40

            if ( I != 0 ) {
               ITMP1 = IWORK( 1 + I-1 )
               W( I ) = W( J )
               IWORK( 1 + I-1 ) = IWORK( 1 + J-1 )
               W( J ) = TMP1
               IWORK( 1 + J-1 ) = ITMP1
               dswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
               if ( INFO != 0 ) {
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               }
            }
         } // 50
      }

      RETURN

      // End of DSBGVX

      }
