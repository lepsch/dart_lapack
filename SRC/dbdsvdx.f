      SUBROUTINE DBDSVDX( UPLO, JOBZ, RANGE, N, D, E, VL, VU, IL, IU, NS, S, Z, LDZ, WORK, IWORK, INFO);

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDZ, N, NS;
      double             VL, VU;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             D( * ), E( * ), S( * ), WORK( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN, HNDRD, MEIGTH;
      const              ZERO = 0.0, ONE = 1.0, TEN = 10.0, HNDRD = 100.0, MEIGTH = -0.1250 ;
      double             FUDGE;
      const              FUDGE = 2.0 ;
      // ..
      // .. Local Scalars ..
      String             RNGVX;
      bool               ALLSV, INDSV, LOWER, SPLIT, SVEQ0, VALSV, WANTZ;
      int                I, ICOLZ, IDBEG, IDEND, IDTGK, IDPTR, IEPTR, IETGK, IIFAIL, IIWORK, ILTGK, IROWU, IROWV, IROWZ, ISBEG, ISPLT, ITEMP, IUTGK, J, K, NTGK, NRU, NRV, NSL;
      double             ABSTOL, EPS, EMIN, MU, NRMU, NRMV, ORTOL, SMAX, SMIN, SQRT2, THRESH, TOL, ULP, VLTGK, VUTGK, ZJTJI;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DDOT, DLAMCH, DNRM2;
      // EXTERNAL IDAMAX, LSAME, DAXPY, DDOT, DLAMCH, DNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSTEVX, DCOPY, DLASET, DSCAL, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      ALLSV = LSAME( RANGE, 'A' );
      VALSV = LSAME( RANGE, 'V' );
      INDSV = LSAME( RANGE, 'I' );
      WANTZ = LSAME( JOBZ, 'V' );
      LOWER = LSAME( UPLO, 'L' );

      INFO = 0;
      if ( !LSAME( UPLO, 'U' ) && !LOWER ) {
         INFO = -1;
      } else if ( !( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( ALLSV || VALSV || INDSV ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( N > 0 ) {
         if ( VALSV ) {
            if ( VL < ZERO ) {
               INFO = -7;
            } else if ( VU <= VL ) {
               INFO = -8;
            }
         } else if ( INDSV ) {
            if ( IL < 1 || IL > MAX( 1, N ) ) {
               INFO = -9;
            } else if ( IU < MIN( N, IL ) || IU > N ) {
               INFO = -10;
            }
         }
      }
      if ( INFO == 0 ) {
         IF( LDZ < 1 || ( WANTZ && LDZ < N*2 ) ) INFO = -14;
      }

      if ( INFO != 0 ) {
         xerbla('DBDSVDX', -INFO );
         return;
      }

      // Quick return if possible (N <= 1)

      NS = 0;
      if (N == 0) RETURN;

      if ( N == 1 ) {
         if ( ALLSV || INDSV ) {
            NS = 1;
            S( 1 ) = ABS( D( 1 ) );
         } else {
            if ( VL < ABS( D( 1 ) ) && VU >= ABS( D( 1 ) ) ) {
               NS = 1;
               S( 1 ) = ABS( D( 1 ) );
            }
         }
         if ( WANTZ ) {
            Z( 1, 1 ) = SIGN( ONE, D( 1 ) );
            Z( 2, 1 ) = ONE;
         }
         return;
      }

      ABSTOL = 2*DLAMCH( 'Safe Minimum' );
      ULP = DLAMCH( 'Precision' );
      EPS = DLAMCH( 'Epsilon' );
      SQRT2 = SQRT( 2.0 );
      ORTOL = SQRT( ULP );

      // Criterion for splitting is taken from DBDSQR when singular
      // values are computed to relative accuracy TOL. (See J. Demmel and
      // W. Kahan, Accurate singular values of bidiagonal matrices, SIAM
      // J. Sci. and Stat. Comput., 11:873–912, 1990.)

      TOL = MAX( TEN, MIN( HNDRD, EPS**MEIGTH ) )*EPS;

      // Compute approximate maximum, minimum singular values.

      I = IDAMAX( N, D, 1 );
      SMAX = ABS( D( I ) );
      I = IDAMAX( N-1, E, 1 );
      SMAX = MAX( SMAX, ABS( E( I ) ) );

      // Compute threshold for neglecting D's and E's.

      SMIN = ABS( D( 1 ) );
      if ( SMIN != ZERO ) {
         MU = SMIN;
         for (I = 2; I <= N; I++) {
            MU = ABS( D( I ) )*( MU / ( MU+ABS( E( I-1 ) ) ) );
            SMIN = MIN( SMIN, MU );
            if (SMIN == ZERO) EXIT;
         }
      }
      SMIN = SMIN / SQRT( DBLE( N ) );
      THRESH = TOL*SMIN;

      // Check for zeros in D and E (splits), i.e. submatrices.

      for (I = 1; I <= N-1; I++) {
         IF( ABS( D( I ) ) <= THRESH ) D( I ) = ZERO;
         IF( ABS( E( I ) ) <= THRESH ) E( I ) = ZERO;
      }
      IF( ABS( D( N ) ) <= THRESH ) D( N ) = ZERO;

      // Pointers for arrays used by DSTEVX.

      IDTGK = 1;
      IETGK = IDTGK + N*2;
      ITEMP = IETGK + N*2;
      IIFAIL = 1;
      IIWORK = IIFAIL + N*2;

      // Set RNGVX, which corresponds to RANGE for DSTEVX in TGK mode.
      // VL,VU or IL,IU are redefined to conform to implementation a)
      // described in the leading comments.

      ILTGK = 0;
      IUTGK = 0;
      VLTGK = ZERO;
      VUTGK = ZERO;

      if ( ALLSV ) {

         // All singular values will be found. We aim at -s (see
         // leading comments) with RNGVX = 'I'. IL and IU are set
         // later (as ILTGK and IUTGK) according to the dimension
         // of the active submatrix.

         RNGVX = 'I';
         if (WANTZ) CALL DLASET( 'F', N*2, N+1, ZERO, ZERO, Z, LDZ );
      } else if ( VALSV ) {

         // Find singular values in a half-open interval. We aim
         // at -s (see leading comments) and we swap VL and VU
         // (as VUTGK and VLTGK), changing their signs.

         RNGVX = 'V';
         VLTGK = -VU;
         VUTGK = -VL;
         WORK( IDTGK:IDTGK+2*N-1 ) = ZERO;
         dcopy(N, D, 1, WORK( IETGK ), 2 );
         dcopy(N-1, E, 1, WORK( IETGK+1 ), 2 );
         dstevx('N', 'V', N*2, WORK( IDTGK ), WORK( IETGK ), VLTGK, VUTGK, ILTGK, ILTGK, ABSTOL, NS, S, Z, LDZ, WORK( ITEMP ), IWORK( IIWORK ), IWORK( IIFAIL ), INFO );
         if ( NS == 0 ) {
            return;
         } else {
            if (WANTZ) CALL DLASET( 'F', N*2, NS, ZERO, ZERO, Z, LDZ );
         }
      } else if ( INDSV ) {

         // Find the IL-th through the IU-th singular values. We aim
         // at -s (see leading comments) and indices are mapped into
         // values, therefore mimicking DSTEBZ, where

         // GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN
         // GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN

         ILTGK = IL;
         IUTGK = IU;
         RNGVX = 'V';
         WORK( IDTGK:IDTGK+2*N-1 ) = ZERO;
         dcopy(N, D, 1, WORK( IETGK ), 2 );
         dcopy(N-1, E, 1, WORK( IETGK+1 ), 2 );
         dstevx('N', 'I', N*2, WORK( IDTGK ), WORK( IETGK ), VLTGK, VLTGK, ILTGK, ILTGK, ABSTOL, NS, S, Z, LDZ, WORK( ITEMP ), IWORK( IIWORK ), IWORK( IIFAIL ), INFO );
         VLTGK = S( 1 ) - FUDGE*SMAX*ULP*N;
         WORK( IDTGK:IDTGK+2*N-1 ) = ZERO;
         dcopy(N, D, 1, WORK( IETGK ), 2 );
         dcopy(N-1, E, 1, WORK( IETGK+1 ), 2 );
         dstevx('N', 'I', N*2, WORK( IDTGK ), WORK( IETGK ), VUTGK, VUTGK, IUTGK, IUTGK, ABSTOL, NS, S, Z, LDZ, WORK( ITEMP ), IWORK( IIWORK ), IWORK( IIFAIL ), INFO );
         VUTGK = S( 1 ) + FUDGE*SMAX*ULP*N;
         VUTGK = MIN( VUTGK, ZERO );

         // If VLTGK=VUTGK, DSTEVX returns an error message,
         // so if needed we change VUTGK slightly.

         if (VLTGK == VUTGK) VLTGK = VLTGK - TOL;

         if (WANTZ) CALL DLASET( 'F', N*2, IU-IL+1, ZERO, ZERO, Z, LDZ);
      }

      // Initialize variables and pointers for S, Z, and WORK.

      // NRU, NRV: number of rows in U and V for the active submatrix
      // IDBEG, ISBEG: offsets for the entries of D and S
      // IROWZ, ICOLZ: offsets for the rows and columns of Z
      // IROWU, IROWV: offsets for the rows of U and V

      NS = 0;
      NRU = 0;
      NRV = 0;
      IDBEG = 1;
      ISBEG = 1;
      IROWZ = 1;
      ICOLZ = 1;
      IROWU = 2;
      IROWV = 1;
      SPLIT = false;
      SVEQ0 = false;

      // Form the tridiagonal TGK matrix.

      S( 1:N ) = ZERO;
      WORK( IETGK+2*N-1 ) = ZERO;
      WORK( IDTGK:IDTGK+2*N-1 ) = ZERO;
      dcopy(N, D, 1, WORK( IETGK ), 2 );
      dcopy(N-1, E, 1, WORK( IETGK+1 ), 2 );


      // Check for splits in two levels, outer level
      // in E and inner level in D.

      DO IEPTR = 2, N*2, 2;
         if ( WORK( IETGK+IEPTR-1 ) == ZERO ) {

            // Split in E (this piece of B is square) or bottom
            // of the (input bidiagonal) matrix.

            ISPLT = IDBEG;
            IDEND = IEPTR - 1;
            DO IDPTR = IDBEG, IDEND, 2;
               if ( WORK( IETGK+IDPTR-1 ) == ZERO ) {

                  // Split in D (rectangular submatrix). Set the number
                  // of rows in U and V (NRU and NRV) accordingly.

                  if ( IDPTR == IDBEG ) {

                     // D=0 at the top.

                     SVEQ0 = true;
                     if ( IDBEG == IDEND) {
                        NRU = 1;
                        NRV = 1;
                     }
                  } else if ( IDPTR == IDEND ) {

                     // D=0 at the bottom.

                     SVEQ0 = true;
                     NRU = (IDEND-ISPLT)/2 + 1;
                     NRV = NRU;
                     if ( ISPLT != IDBEG ) {
                        NRU = NRU + 1;
                     }
                  } else {
                     if ( ISPLT == IDBEG ) {

                        // Split: top rectangular submatrix.

                        NRU = (IDPTR-IDBEG)/2;
                        NRV = NRU + 1;
                     } else {

                        // Split: middle square submatrix.

                        NRU = (IDPTR-ISPLT)/2 + 1;
                        NRV = NRU;
                     }
                  }
               } else if ( IDPTR == IDEND ) {

                  // Last entry of D in the active submatrix.

                  if ( ISPLT == IDBEG ) {

                     // No split (trivial case).

                     NRU = (IDEND-IDBEG)/2 + 1;
                     NRV = NRU;
                  } else {

                     // Split: bottom rectangular submatrix.

                     NRV = (IDEND-ISPLT)/2 + 1;
                     NRU = NRV + 1;
                  }
               }

               NTGK = NRU + NRV;

               if ( NTGK > 0 ) {

                  // Compute eigenvalues/vectors of the active
                  // submatrix according to RANGE:
                  // if RANGE='A' (ALLSV) then RNGVX = 'I'
                  // if RANGE='V' (VALSV) then RNGVX = 'V'
                  // if RANGE='I' (INDSV) then RNGVX = 'V'

                  ILTGK = 1;
                  IUTGK = NTGK / 2;
                  if ( ALLSV || VUTGK == ZERO ) {
                     if ( SVEQ0 || SMIN < EPS || MOD(NTGK,2) > 0 ) {
                         // Special case: eigenvalue equal to zero or very
                         // small, additional eigenvector is needed.
                         IUTGK = IUTGK + 1;
                     }
                  }

                  // Workspace needed by DSTEVX:
                  // WORK( ITEMP: ): 2*5*NTGK
                  // IWORK( 1: ): 2*6*NTGK

                  dstevx(JOBZ, RNGVX, NTGK, WORK( IDTGK+ISPLT-1 ), WORK( IETGK+ISPLT-1 ), VLTGK, VUTGK, ILTGK, IUTGK, ABSTOL, NSL, S( ISBEG ), Z( IROWZ,ICOLZ ), LDZ, WORK( ITEMP ), IWORK( IIWORK ), IWORK( IIFAIL ), INFO );
                  if ( INFO != 0 ) {
                     // Exit with the error code from DSTEVX.
                     return;
                  }
                  EMIN = ABS( MAXVAL( S( ISBEG:ISBEG+NSL-1 ) ) );

                  if ( NSL > 0 && WANTZ ) {

                     // Normalize u=Z([2,4,...],:) and v=Z([1,3,...],:),
                     // changing the sign of v as discussed in the leading
                     // comments. The norms of u and v may be (slightly)
                     // different from 1/sqrt(2) if the corresponding
                     // eigenvalues are very small or too close. We check
                     // those norms and, if needed, reorthogonalize the
                     // vectors.

                     if ( NSL > 1 && VUTGK == ZERO && MOD(NTGK,2) == 0 && EMIN == 0 && !SPLIT ) {

                        // D=0 at the top or bottom of the active submatrix:
                        // one eigenvalue is equal to zero; concatenate the
                        // eigenvectors corresponding to the two smallest
                        // eigenvalues.

                        Z( IROWZ:IROWZ+NTGK-1,ICOLZ+NSL-2 ) = Z( IROWZ:IROWZ+NTGK-1,ICOLZ+NSL-2 ) + Z( IROWZ:IROWZ+NTGK-1,ICOLZ+NSL-1 );
                        Z( IROWZ:IROWZ+NTGK-1,ICOLZ+NSL-1 ) = ZERO;
                        // IF( IUTGK*2 > NTGK ) THEN
                           // Eigenvalue equal to zero or very small.
                           // NSL = NSL - 1
                        // END IF
                     }

                     DO I = 0, MIN( NSL-1, NRU-1 );
                        NRMU = DNRM2( NRU, Z( IROWU, ICOLZ+I ), 2 );
                        if ( NRMU == ZERO ) {
                           INFO = N*2 + 1;
                           return;
                        }
                        CALL DSCAL( NRU, ONE/NRMU, Z( IROWU,ICOLZ+I ), 2 )                         IF( NRMU != ONE && ABS( NRMU-ORTOL )*SQRT2 > ONE ) THEN;
                           for (J = 0; J <= I-1; J++) {
                              ZJTJI = -DDOT( NRU, Z( IROWU, ICOLZ+J ), 2, Z( IROWU, ICOLZ+I ), 2 );
                              daxpy(NRU, ZJTJI, Z( IROWU, ICOLZ+J ), 2, Z( IROWU, ICOLZ+I ), 2 );
                           }
                           NRMU = DNRM2( NRU, Z( IROWU, ICOLZ+I ), 2 );
                           dscal(NRU, ONE/NRMU, Z( IROWU,ICOLZ+I ), 2 );
                        }
                     }
                     DO I = 0, MIN( NSL-1, NRV-1 );
                        NRMV = DNRM2( NRV, Z( IROWV, ICOLZ+I ), 2 );
                        if ( NRMV == ZERO ) {
                           INFO = N*2 + 1;
                           return;
                        }
                        CALL DSCAL( NRV, -ONE/NRMV, Z( IROWV,ICOLZ+I ), 2 )                         IF( NRMV != ONE && ABS( NRMV-ORTOL )*SQRT2 > ONE ) THEN;
                           for (J = 0; J <= I-1; J++) {
                              ZJTJI = -DDOT( NRV, Z( IROWV, ICOLZ+J ), 2, Z( IROWV, ICOLZ+I ), 2 );
                              daxpy(NRU, ZJTJI, Z( IROWV, ICOLZ+J ), 2, Z( IROWV, ICOLZ+I ), 2 );
                           }
                           NRMV = DNRM2( NRV, Z( IROWV, ICOLZ+I ), 2 );
                           dscal(NRV, ONE/NRMV, Z( IROWV,ICOLZ+I ), 2 );
                        }
                     }
                     if ( VUTGK == ZERO && IDPTR < IDEND && MOD(NTGK,2) > 0 ) {

                        // D=0 in the middle of the active submatrix (one
                        // eigenvalue is equal to zero): save the corresponding
                        // eigenvector for later use (when bottom of the
                        // active submatrix is reached).

                        SPLIT = true;
                        Z( IROWZ:IROWZ+NTGK-1,N+1 ) = Z( IROWZ:IROWZ+NTGK-1,NS+NSL )                         Z( IROWZ:IROWZ+NTGK-1,NS+NSL ) = ZERO;
                     }
                  END IF !** WANTZ **!;

                  NSL = MIN( NSL, NRU );
                  SVEQ0 = false;

                  // Absolute values of the eigenvalues of TGK.

                  for (I = 0; I <= NSL-1; I++) {
                     S( ISBEG+I ) = ABS( S( ISBEG+I ) );
                  }

                  // Update pointers for TGK, S and Z.

                  ISBEG = ISBEG + NSL;
                  IROWZ = IROWZ + NTGK;
                  ICOLZ = ICOLZ + NSL;
                  IROWU = IROWZ;
                  IROWV = IROWZ + 1;
                  ISPLT = IDPTR + 1;
                  NS = NS + NSL;
                  NRU = 0;
                  NRV = 0;
               END IF !** NTGK > 0 **!;
               if ( IROWZ < N*2 && WANTZ ) {
                  Z( 1:IROWZ-1, ICOLZ ) = ZERO;
               }
            END DO !** IDPTR loop **!;
            if ( SPLIT && WANTZ ) {

               // Bring back eigenvector corresponding
               // to eigenvalue equal to zero.

               Z( IDBEG:IDEND-NTGK+1,ISBEG-1 ) = Z( IDBEG:IDEND-NTGK+1,ISBEG-1 ) + Z( IDBEG:IDEND-NTGK+1,N+1 );
               Z( IDBEG:IDEND-NTGK+1,N+1 ) = 0;
            }
            IROWV = IROWV - 1;
            IROWU = IROWU + 1;
            IDBEG = IEPTR + 1;
            SVEQ0 = false;
            SPLIT = false;
         END IF !** Check for split in E **!;
      END DO !** IEPTR loop **!;

      // Sort the singular values into decreasing order (insertion sort on
      // singular values, but only one transposition per singular vector)

      for (I = 1; I <= NS-1; I++) {
         K = 1;
         SMIN = S( 1 );
         for (J = 2; J <= NS + 1 - I; J++) {
            if ( S( J ) <= SMIN ) {
               K = J;
               SMIN = S( J );
            }
         }
         if ( K != NS+1-I ) {
            S( K ) = S( NS+1-I );
            S( NS+1-I ) = SMIN;
            if (WANTZ) CALL DSWAP( N*2, Z( 1,K ), 1, Z( 1,NS+1-I ), 1 );
         }
      }

      // If RANGE=I, check for singular values/vectors to be discarded.

      if ( INDSV ) {
         K = IU - IL + 1;
         if ( K < NS ) {
            S( K+1:NS ) = ZERO;
            if (WANTZ) Z( 1:N*2,K+1:NS ) = ZERO;
            NS = K;
         }
      }

      // Reorder Z: U = Z( 1:N,1:NS ), V = Z( N+1:N*2,1:NS ).
      // If B is a lower diagonal, swap U and V.

      if ( WANTZ ) {
      for (I = 1; I <= NS; I++) {
         dcopy(N*2, Z( 1,I ), 1, WORK, 1 );
         if ( LOWER ) {
            dcopy(N, WORK( 2 ), 2, Z( N+1,I ), 1 );
            dcopy(N, WORK( 1 ), 2, Z( 1  ,I ), 1 );
         } else {
            dcopy(N, WORK( 2 ), 2, Z( 1  ,I ), 1 );
            dcopy(N, WORK( 1 ), 2, Z( N+1,I ), 1 );
         }
      }
      }

      return;

      // End of DBDSVDX

      }
