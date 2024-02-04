      void zlatrs3(UPLO, TRANS, DIAG, NORMIN, N, NRHS, A, LDA, X, LDX, SCALE, CNORM, WORK, LWORK, INFO ) {
      // IMPLICIT NONE

      // .. Scalar Arguments ..
      String             DIAG, TRANS, NORMIN, UPLO;
      int                INFO, LDA, LWORK, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      Complex         A( LDA, * ), X( LDX, * );
      double             CNORM( * ), SCALE( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      const              CZERO = ( 0.0, 0.0 ) ;
      int                NBMAX, NBMIN, NBRHS, NRHSMIN;
      const              NRHSMIN = 2, NBRHS = 32 ;
      const              NBMIN = 8, NBMAX = 64 ;
      // ..
      // .. Local Arrays ..
      double             W( NBMAX ), XNRM( NBRHS );
      // ..
      // .. Local Scalars ..
      bool               LQUERY, NOTRAN, NOUNIT, UPPER;
      int                AWRK, I, IFIRST, IINC, ILAST, II, I1, I2, J, JFIRST, JINC, JLAST, J1, J2, K, KK, K1, K2, LANRM, LDS, LSCALE, NB, NBA, NBX, RHS, LWMIN;
      double             ANRM, BIGNUM, BNRM, RSCAL, SCAL, SCALOC, SCAMIN, SMLNUM, TMAX;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- double             DLAMCH, ZLANGE, DLARMM;
      // EXTERNAL ILAENV, lsame, DLAMCH, ZLANGE, DLARMM
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLATRS, ZDSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      NOTRAN = lsame( TRANS, 'N' );
      NOUNIT = lsame( DIAG, 'N' );
      LQUERY = ( LWORK == -1 );

      // Partition A and X into blocks.

      NB = max( NBMIN, ILAENV( 1, 'ZLATRS', '', N, N, -1, -1 ) );
      NB = min( NBMAX, NB );
      NBA = max( 1, (N + NB - 1) / NB );
      NBX = max( 1, (NRHS + NBRHS - 1) / NBRHS );

      // Compute the workspace

      // The workspace comprises two parts.
      // The first part stores the local scale factors. Each simultaneously
      // computed right-hand side requires one local scale factor per block
      // row. WORK( I + KK * LDS ) is the scale factor of the vector
      // segment associated with the I-th block row and the KK-th vector
      // in the block column.

      LSCALE = NBA * max( NBA, min( NRHS, NBRHS ) );
      LDS = NBA;

      // The second part stores upper bounds of the triangular A. There are
      // a total of NBA x NBA blocks, of which only the upper triangular
      // part or the lower triangular part is referenced. The upper bound of
      // the block A( I, J ) is stored as WORK( AWRK + I + J * NBA ).

      LANRM = NBA * NBA;
      AWRK = LSCALE;

      if ( min( N, NRHS ) == 0 ) {
         LWMIN = 1;
      } else {
         LWMIN = LSCALE + LANRM;
      }
      WORK[1] = LWMIN;

      // Test the input parameters.

      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !NOTRAN && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( !NOUNIT && !lsame( DIAG, 'U' ) ) {
         INFO = -3;
      } else if ( !lsame( NORMIN, 'Y' ) && !lsame( NORMIN, 'N' ) ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( NRHS < 0 ) {
         INFO = -6;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -8;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -10;
      } else if ( !LQUERY && LWORK < LWMIN ) {
         INFO = -14;
      }
      if ( INFO != 0 ) {
         xerbla('ZLATRS3', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Initialize scaling factors

      for (KK = 1; KK <= NRHS; KK++) {
         SCALE[KK] = ONE;
      }

      // Quick return if possible

      if( min( N, NRHS ) == 0 ) return;

      // Determine machine dependent constant to control overflow.

      BIGNUM = DLAMCH( 'Overflow' );
      SMLNUM = DLAMCH( 'Safe Minimum' );

      // Use unblocked code for small problems

      if ( NRHS < NRHSMIN ) {
         zlatrs(UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X( 1, 1), SCALE( 1 ), CNORM, INFO );
         for (K = 2; K <= NRHS; K++) {
            zlatrs(UPLO, TRANS, DIAG, 'Y', N, A, LDA, X( 1, K ), SCALE( K ), CNORM, INFO );
         }
         return;
      }

      // Compute norms of blocks of A excluding diagonal blocks and find
      // the block with the largest norm TMAX.

      TMAX = ZERO;
      for (J = 1; J <= NBA; J++) {
         J1 = (J-1)*NB + 1;
         J2 = min( J*NB, N ) + 1;
         if ( UPPER ) {
            IFIRST = 1;
            ILAST = J - 1;
         } else {
            IFIRST = J + 1;
            ILAST = NBA;
         }
         for (I = IFIRST; I <= ILAST; I++) {
            I1 = (I-1)*NB + 1;
            I2 = min( I*NB, N ) + 1;

            // Compute upper bound of A( I1:I2-1, J1:J2-1 ).

            if ( NOTRAN ) {
               ANRM = ZLANGE( 'I', I2-I1, J2-J1, A( I1, J1 ), LDA, W );
               WORK[AWRK + I+(J-1)*NBA] = ANRM;
            } else {
               ANRM = ZLANGE( '1', I2-I1, J2-J1, A( I1, J1 ), LDA, W );
               WORK[AWRK + J+(I-1) * NBA] = ANRM;
            }
            TMAX = max( TMAX, ANRM );
         }
      }

      if ( !TMAX <= DLAMCH('Overflow') ) {

         // Some matrix entries have huge absolute value. At least one upper
         // bound norm( A(I1:I2-1, J1:J2-1), 'I') is not a valid floating-point
         // number, either due to overflow in LANGE or due to Inf in A.
         // Fall back to LATRS. Set normin = 'N' for every right-hand side to
         // force computation of TSCAL in LATRS to avoid the likely overflow
         // in the computation of the column norms CNORM.

         for (K = 1; K <= NRHS; K++) {
            zlatrs(UPLO, TRANS, DIAG, 'N', N, A, LDA, X( 1, K ), SCALE( K ), CNORM, INFO );
         }
         return;
      }

      // Every right-hand side requires workspace to store NBA local scale
      // factors. To save workspace, X is computed successively in block columns
      // of width NBRHS, requiring a total of NBA x NBRHS space. If sufficient
      // workspace is available, larger values of NBRHS or NBRHS = NRHS are viable.
      for (K = 1; K <= NBX; K++) {
         // Loop over block columns (index = K) of X and, for column-wise scalings,
         // over individual columns (index = KK).
         // K1: column index of the first column in X( J, K )
         // K2: column index of the first column in X( J, K+1 )
         // so the K2 - K1 is the column count of the block X( J, K )
         K1 = (K-1)*NBRHS + 1;
         K2 = min( K*NBRHS, NRHS ) + 1;

         // Initialize local scaling factors of current block column X( J, K )

         for (KK = 1; KK <= K2 - K1; KK++) {
            for (I = 1; I <= NBA; I++) {
               WORK[I+KK*LDS] = ONE;
            }
         }

         if ( NOTRAN ) {

            // Solve A * X(:, K1:K2-1) = B * diag(scale(K1:K2-1))

            if ( UPPER ) {
               JFIRST = NBA;
               JLAST = 1;
               JINC = -1;
            } else {
               JFIRST = 1;
               JLAST = NBA;
               JINC = 1;
            }
         } else {

            // Solve op(A) * X(:, K1:K2-1) = B * diag(scale(K1:K2-1))
            // where op(A) = A**T or op(A) = A**H

            if ( UPPER ) {
               JFIRST = 1;
               JLAST = NBA;
               JINC = 1;
            } else {
               JFIRST = NBA;
               JLAST = 1;
               JINC = -1;
            }
         }

         for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
            // J1: row index of the first row in A( J, J )
            // J2: row index of the first row in A( J+1, J+1 )
            // so that J2 - J1 is the row count of the block A( J, J )
            J1 = (J-1)*NB + 1;
            J2 = min( J*NB, N ) + 1;

            // Solve op(A( J, J )) * X( J, RHS ) = SCALOC * B( J, RHS )

            for (KK = 1; KK <= K2 - K1; KK++) {
               RHS = K1 + KK - 1;
               if ( KK == 1 ) {
                  zlatrs(UPLO, TRANS, DIAG, 'N', J2-J1, A( J1, J1 ), LDA, X( J1, RHS ), SCALOC, CNORM, INFO );
               } else {
                  zlatrs(UPLO, TRANS, DIAG, 'Y', J2-J1, A( J1, J1 ), LDA, X( J1, RHS ), SCALOC, CNORM, INFO );
               }
               // Find largest absolute value entry in the vector segment
               // X( J1:J2-1, RHS ) as an upper bound for the worst case
               // growth in the linear updates.
               XNRM[KK] = ZLANGE( 'I', J2-J1, 1, X( J1, RHS ), LDX, W );

               if ( SCALOC == ZERO ) {
                  // LATRS found that A is singular through A(j,j) = 0.
                  // Reset the computation x(1:n) = 0, x(j) = 1, SCALE = 0
                  // and compute op(A)*x = 0. Note that X(J1:J2-1, KK) is
                  // set by LATRS.
                  SCALE[RHS] = ZERO;
                  for (II = 1; II <= J1-1; II++) {
                     X[II, KK] = CZERO;
                  }
                  for (II = J2; II <= N; II++) {
                     X[II, KK] = CZERO;
                  }
                  // Discard the local scale factors.
                  for (II = 1; II <= NBA; II++) {
                     WORK[II+KK*LDS] = ONE;
                  }
                  SCALOC = ONE;
               } else if ( SCALOC*WORK( J+KK*LDS ) == ZERO ) {
                  // LATRS computed a valid scale factor, but combined with
                  // the current scaling the solution does not have a
                  // scale factor > 0.

                  // Set WORK( J+KK*LDS ) to smallest valid scale
                  // factor and increase SCALOC accordingly.
                  SCAL = WORK( J+KK*LDS ) / SMLNUM;
                  SCALOC = SCALOC * SCAL;
                  WORK[J+KK*LDS] = SMLNUM;
                  // If LATRS overestimated the growth, x may be
                  // rescaled to preserve a valid combined scale
                  // factor WORK( J, KK ) > 0.
                  RSCAL = ONE / SCALOC;
                  if ( XNRM( KK )*RSCAL <= BIGNUM ) {
                     XNRM[KK] = XNRM( KK ) * RSCAL;
                     zdscal(J2-J1, RSCAL, X( J1, RHS ), 1 );
                     SCALOC = ONE;
                  } else {
                     // The system op(A) * x = b is badly scaled and its
                     // solution cannot be represented as (1/scale) * x.
                     // Set x to zero. This approach deviates from LATRS
                     // where a completely meaningless non-zero vector
                     // is returned that is not a solution to op(A) * x = b.
                     SCALE[RHS] = ZERO;
                     for (II = 1; II <= N; II++) {
                        X[II, KK] = CZERO;
                     }
                     // Discard the local scale factors.
                     for (II = 1; II <= NBA; II++) {
                        WORK[II+KK*LDS] = ONE;
                     }
                     SCALOC = ONE;
                  }
               }
               SCALOC = SCALOC * WORK( J+KK*LDS );
               WORK[J+KK*LDS] = SCALOC;
            }

            // Linear block updates

            if ( NOTRAN ) {
               if ( UPPER ) {
                  IFIRST = J - 1;
                  ILAST = 1;
                  IINC = -1;
               } else {
                  IFIRST = J + 1;
                  ILAST = NBA;
                  IINC = 1;
               }
            } else {
               if ( UPPER ) {
                  IFIRST = J + 1;
                  ILAST = NBA;
                  IINC = 1;
               } else {
                  IFIRST = J - 1;
                  ILAST = 1;
                  IINC = -1;
               }
            }

            for (I = IFIRST; IINC < 0 ? I >= ILAST : I <= ILAST; I += IINC) {
               // I1: row index of the first column in X( I, K )
               // I2: row index of the first column in X( I+1, K )
               // so the I2 - I1 is the row count of the block X( I, K )
               I1 = (I-1)*NB + 1;
               I2 = min( I*NB, N ) + 1;

               // Prepare the linear update to be executed with GEMM.
               // For each column, compute a consistent scaling, a
               // scaling factor to survive the linear update, and
               // rescale the column segments, if necessary. Then
               // the linear update is safely executed.

               for (KK = 1; KK <= K2 - K1; KK++) {
                  RHS = K1 + KK - 1;
                  // Compute consistent scaling
                  SCAMIN = min( WORK( I+KK*LDS), WORK( J+KK*LDS ) );

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  BNRM = ZLANGE( 'I', I2-I1, 1, X( I1, RHS ), LDX, W );
                  BNRM = BNRM*( SCAMIN / WORK( I+KK*LDS ) );
                  XNRM[KK] = XNRM( KK )*( SCAMIN / WORK( J+KK*LDS) );
                  ANRM = WORK( AWRK + I+(J-1)*NBA );
                  SCALOC = DLARMM( ANRM, XNRM( KK ), BNRM );

                  // Simultaneously apply the robust update factor and the
                  // consistency scaling factor to X( I, KK ) and X( J, KK ).

                  SCAL = ( SCAMIN / WORK( I+KK*LDS) )*SCALOC;
                  if ( SCAL != ONE ) {
                     zdscal(I2-I1, SCAL, X( I1, RHS ), 1 );
                     WORK[I+KK*LDS] = SCAMIN*SCALOC;
                  }

                  SCAL = ( SCAMIN / WORK( J+KK*LDS ) )*SCALOC;
                  if ( SCAL != ONE ) {
                     zdscal(J2-J1, SCAL, X( J1, RHS ), 1 );
                     WORK[J+KK*LDS] = SCAMIN*SCALOC;
                  }
               }

               if ( NOTRAN ) {

                  // B( I, K ) := B( I, K ) - A( I, J ) * X( J, K )

                  zgemm('N', 'N', I2-I1, K2-K1, J2-J1, -CONE, A( I1, J1 ), LDA, X( J1, K1 ), LDX, CONE, X( I1, K1 ), LDX );
               } else if ( lsame( TRANS, 'T' ) ) {

                  // B( I, K ) := B( I, K ) - A( I, J )**T * X( J, K )

                  zgemm('T', 'N', I2-I1, K2-K1, J2-J1, -CONE, A( J1, I1 ), LDA, X( J1, K1 ), LDX, CONE, X( I1, K1 ), LDX );
               } else {

                  // B( I, K ) := B( I, K ) - A( I, J )**H * X( J, K )

                  zgemm('C', 'N', I2-I1, K2-K1, J2-J1, -CONE, A( J1, I1 ), LDA, X( J1, K1 ), LDX, CONE, X( I1, K1 ), LDX );
               }
            }
         }


         // Reduce local scaling factors

         for (KK = 1; KK <= K2 - K1; KK++) {
            RHS = K1 + KK - 1;
            for (I = 1; I <= NBA; I++) {
               SCALE[RHS] = min( SCALE( RHS ), WORK( I+KK*LDS ) );
            }
         }

         // Realize consistent scaling

         for (KK = 1; KK <= K2 - K1; KK++) {
            RHS = K1 + KK - 1;
            if ( SCALE( RHS ) != ONE && SCALE( RHS ) != ZERO ) {
               for (I = 1; I <= NBA; I++) {
                  I1 = (I - 1) * NB + 1;
                  I2 = min( I * NB, N ) + 1;
                  SCAL = SCALE( RHS ) / WORK( I+KK*LDS );
                  if (SCAL != ONE) zdscal( I2-I1, SCAL, X( I1, RHS ), 1 );
               }
            }
         }
      }
      return;
      }
