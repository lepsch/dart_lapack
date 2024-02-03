      SUBROUTINE ZTRSYL3( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, LDC, SCALE, SWORK, LDSWORK, INFO )
      IMPLICIT NONE

      // .. Scalar Arguments ..
      String             TRANA, TRANB;
      int                INFO, ISGN, LDA, LDB, LDC, LDSWORK, M, N;
      double             SCALE;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
      double             SWORK( LDSWORK, * );
      // ..
      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      COMPLEX*16         CONE
      const              CONE = ( 1.0D0, 0.0D0 ) ;
      // ..
      // .. Local Scalars ..
      bool               NOTRNA, NOTRNB, LQUERY;
      int                AWRK, BWRK, I, I1, I2, IINFO, J, J1, J2, JJ, K, K1, K2, L, L1, L2, LL, NBA, NB, NBB;
      double             ANRM, BIGNUM, BNRM, CNRM, SCAL, SCALOC, SCAMIN, SGN, XNRM, BUF, SMLNUM;
      COMPLEX*16         CSGN
      // ..
      // .. Local Arrays ..
      double             WNRM( MAX( M, N ) );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, DLARMM, ZLANGE;
      // EXTERNAL DLAMCH, DLARMM, ILAENV, LSAME, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDSCAL, ZGEMM, ZLASCL, ZTRSYL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, EXPONENT, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Decode and Test input parameters

      NOTRNA = LSAME( TRANA, 'N' )
      NOTRNB = LSAME( TRANB, 'N' )

      // Use the same block size for all matrices.

      NB = MAX( 8, ILAENV( 1, 'ZTRSYL', '', M, N, -1, -1) )

      // Compute number of blocks in A and B

      NBA = MAX( 1, (M + NB - 1) / NB )
      NBB = MAX( 1, (N + NB - 1) / NB )

      // Compute workspace

      INFO = 0
      LQUERY = ( LDSWORK.EQ.-1 )
      if ( LQUERY ) {
         LDSWORK = 2
         SWORK(1,1) = MAX( NBA, NBB )
         SWORK(2,1) = 2 * NBB + NBA
      }

      // Test the input arguments

      if ( .NOT.NOTRNA .AND. .NOT. LSAME( TRANA, 'C' ) ) {
         INFO = -1
      } else if ( .NOT.NOTRNB .AND. .NOT. LSAME( TRANB, 'C' ) ) {
         INFO = -2
      } else if ( ISGN.NE.1 .AND. ISGN.NE.-1 ) {
         INFO = -3
      } else if ( M.LT.0 ) {
         INFO = -4
      } else if ( N.LT.0 ) {
         INFO = -5
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -7
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -9
      } else if ( LDC.LT.MAX( 1, M ) ) {
         INFO = -11
      }
      if ( INFO.NE.0 ) {
         xerbla('ZTRSYL3', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      SCALE = ONE
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      // Use unblocked code for small problems or if insufficient
      // workspace is provided

      if ( MIN( NBA, NBB ).EQ.1 .OR. LDSWORK.LT.MAX( NBA, NBB ) ) {
        ztrsyl(TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, LDC, SCALE, INFO );
        RETURN
      }

      // Set constants to control overflow

      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM

      // Set local scaling factors.

      for (L = 1; L <= NBB; L++) {
         for (K = 1; K <= NBA; K++) {
            SWORK( K, L ) = ONE
         END DO
      END DO

      // Fallback scaling factor to prevent flushing of SWORK( K, L ) to zero.
      // This scaling is to ensure compatibility with TRSYL and may get flushed.

      BUF = ONE

       // Compute upper bounds of blocks of A and B

      AWRK = NBB
      for (K = 1; K <= NBA; K++) {
         K1 = (K - 1) * NB + 1
         K2 = MIN( K * NB, M ) + 1
         for (L = K; L <= NBA; L++) {
            L1 = (L - 1) * NB + 1
            L2 = MIN( L * NB, M ) + 1
            if ( NOTRNA ) {
               SWORK( K, AWRK + L ) = ZLANGE( 'I', K2-K1, L2-L1, A( K1, L1 ), LDA, WNRM )
            } else {
               SWORK( L, AWRK + K ) = ZLANGE( '1', K2-K1, L2-L1, A( K1, L1 ), LDA, WNRM )
            }
         END DO
      END DO
      BWRK = NBB + NBA
      for (K = 1; K <= NBB; K++) {
         K1 = (K - 1) * NB + 1
         K2 = MIN( K * NB, N ) + 1
         for (L = K; L <= NBB; L++) {
            L1 = (L - 1) * NB + 1
            L2 = MIN( L * NB, N ) + 1
            if ( NOTRNB ) {
               SWORK( K, BWRK + L ) = ZLANGE( 'I', K2-K1, L2-L1, B( K1, L1 ), LDB, WNRM )
            } else {
               SWORK( L, BWRK + K ) = ZLANGE( '1', K2-K1, L2-L1, B( K1, L1 ), LDB, WNRM )
            }
         END DO
      END DO

      SGN = DBLE( ISGN )
      CSGN = DCMPLX( SGN, ZERO )

      if ( NOTRNA .AND. NOTRNB ) {

         // Solve    A*X + ISGN*X*B = scale*C.

         // The (K,L)th block of X is determined starting from
         // bottom-left corner column by column by

          // A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

         // Where
                   // M                         L-1
         // R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
                 // I=K+1                       J=1

         // Start loop over block rows (index = K) and block columns (index = L)

         DO K = NBA, 1, -1

            // K1: row index of the first row in X( K, L )
            // K2: row index of the first row in X( K+1, L )
            // so the K2 - K1 is the column count of the block X( K, L )

            K1 = (K - 1) * NB + 1
            K2 = MIN( K * NB, M ) + 1
            for (L = 1; L <= NBB; L++) {

               // L1: column index of the first column in X( K, L )
               // L2: column index of the first column in X( K, L + 1)
               // so that L2 - L1 is the row count of the block X( K, L )

               L1 = (L - 1) * NB + 1
               L2 = MIN( L * NB, N ) + 1

               ztrsyl(TRANA, TRANB, ISGN, K2-K1, L2-L1, A( K1, K1 ), LDA, B( L1, L1 ), LDB, C( K1, L1 ), LDC, SCALOC, IINFO );
               INFO = MAX( INFO, IINFO )

               if ( SCALOC * SWORK( K, L ) .EQ. ZERO ) {
                  if ( SCALOC .EQ. ZERO ) {
                     // The magnitude of the largest entry of X(K1:K2-1, L1:L2-1)
                     // is larger than the product of BIGNUM**2 and cannot be
                     // represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1).
                     // Mark the computation as pointless.
                     BUF = ZERO
                  } else {
                     BUF = BUF*2.D0**EXPONENT( SCALOC )
                  }
                  for (JJ = 1; JJ <= NBB; JJ++) {
                     for (LL = 1; LL <= NBA; LL++) {
                        // Bound by BIGNUM to not introduce Inf. The value
                        // is irrelevant; corresponding entries of the
                        // solution will be flushed in consistency scaling.
                        SWORK( LL, JJ ) = MIN( BIGNUM, SWORK( LL, JJ ) / 2.D0**EXPONENT( SCALOC ) )
                     END DO
                  END DO
               }
               SWORK( K, L ) = SCALOC * SWORK( K, L )
               XNRM = ZLANGE( 'I', K2-K1, L2-L1, C( K1, L1 ), LDC, WNRM )

               DO I = K - 1, 1, -1

                  // C( I, L ) := C( I, L ) - A( I, K ) * C( K, L )

                  I1 = (I - 1) * NB + 1
                  I2 = MIN( I * NB, M ) + 1

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = ZLANGE( 'I', I2-I1, L2-L1, C( I1, L1 ), LDC, WNRM )
                  SCAMIN = MIN( SWORK( I, L ), SWORK( K, L ) )
                  CNRM = CNRM * ( SCAMIN / SWORK( I, L ) )
                  XNRM = XNRM * ( SCAMIN / SWORK( K, L ) )
                  ANRM = SWORK( I, AWRK + K )
                  SCALOC = DLARMM( ANRM, XNRM, CNRM )
                  if ( SCALOC * SCAMIN .EQ. ZERO ) {
                     // Use second scaling factor to prevent flushing to zero.
                     BUF = BUF*2.D0**EXPONENT( SCALOC )
                     for (JJ = 1; JJ <= NBB; JJ++) {
                        for (LL = 1; LL <= NBA; LL++) {
                        SWORK( LL, JJ ) = MIN( BIGNUM, SWORK( LL, JJ ) / 2.D0**EXPONENT( SCALOC ) )
                        END DO
                     END DO
                     SCAMIN = SCAMIN / 2.D0**EXPONENT( SCALOC )
                     SCALOC = SCALOC / 2.D0**EXPONENT( SCALOC )
                  }
                  CNRM = CNRM * SCALOC
                  XNRM = XNRM * SCALOC

                  // Simultaneously apply the robust update factor and the
                  // consistency scaling factor to C( I, L ) and C( K, L ).

                  SCAL = ( SCAMIN / SWORK( K, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                      DO JJ = L1, L2-1
                         zdscal(K2-K1, SCAL, C( K1, JJ ), 1);
                      END DO
                  ENDIF

                  SCAL = ( SCAMIN / SWORK( I, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                      DO LL = L1, L2-1
                         zdscal(I2-I1, SCAL, C( I1, LL ), 1);
                      END DO
                  ENDIF

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( I, L ) = SCAMIN * SCALOC

                  zgemm('N', 'N', I2-I1, L2-L1, K2-K1, -CONE, A( I1, K1 ), LDA, C( K1, L1 ), LDC, CONE, C( I1, L1 ), LDC );

               END DO

               DO J = L + 1, NBB

                  // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( L, J )

                  J1 = (J - 1) * NB + 1
                  J2 = MIN( J * NB, N ) + 1

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = ZLANGE( 'I', K2-K1, J2-J1, C( K1, J1 ), LDC, WNRM )
                  SCAMIN = MIN( SWORK( K, J ), SWORK( K, L ) )
                  CNRM = CNRM * ( SCAMIN / SWORK( K, J ) )
                  XNRM = XNRM * ( SCAMIN / SWORK( K, L ) )
                  BNRM = SWORK(L, BWRK + J)
                  SCALOC = DLARMM( BNRM, XNRM, CNRM )
                  if ( SCALOC * SCAMIN .EQ. ZERO ) {
                     // Use second scaling factor to prevent flushing to zero.
                     BUF = BUF*2.D0**EXPONENT( SCALOC )
                     for (JJ = 1; JJ <= NBB; JJ++) {
                        for (LL = 1; LL <= NBA; LL++) {
                        SWORK( LL, JJ ) = MIN( BIGNUM, SWORK( LL, JJ ) / 2.D0**EXPONENT( SCALOC ) )
                        END DO
                     END DO
                     SCAMIN = SCAMIN / 2.D0**EXPONENT( SCALOC )
                     SCALOC = SCALOC / 2.D0**EXPONENT( SCALOC )
                  }
                  CNRM = CNRM * SCALOC
                  XNRM = XNRM * SCALOC

                  // Simultaneously apply the robust update factor and the
                  // consistency scaling factor to C( K, J ) and C( K, L).

                  SCAL = ( SCAMIN / SWORK( K, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     DO LL = L1, L2-1
                        zdscal(K2-K1, SCAL, C( K1, LL ), 1 );
                     END DO
                  ENDIF

                  SCAL = ( SCAMIN / SWORK( K, J ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                      DO JJ = J1, J2-1
                         zdscal(K2-K1, SCAL, C( K1, JJ ), 1 );
                      END DO
                  ENDIF

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( K, J ) = SCAMIN * SCALOC

                  zgemm('N', 'N', K2-K1, J2-J1, L2-L1, -CSGN, C( K1, L1 ), LDC, B( L1, J1 ), LDB, CONE, C( K1, J1 ), LDC );
               END DO
            END DO
         END DO
      } else if ( .NOT.NOTRNA .AND. NOTRNB ) {

         // Solve    A**H *X + ISGN*X*B = scale*C.

         // The (K,L)th block of X is determined starting from
         // upper-left corner column by column by

           // A(K,K)**H*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

         // Where
                    // K-1                        L-1
           // R(K,L) = SUM [A(I,K)**H*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
                    // I=1                        J=1

         // Start loop over block rows (index = K) and block columns (index = L)

         for (K = 1; K <= NBA; K++) {

            // K1: row index of the first row in X( K, L )
            // K2: row index of the first row in X( K+1, L )
            // so the K2 - K1 is the column count of the block X( K, L )

            K1 = (K - 1) * NB + 1
            K2 = MIN( K * NB, M ) + 1
            for (L = 1; L <= NBB; L++) {

               // L1: column index of the first column in X( K, L )
               // L2: column index of the first column in X( K, L + 1)
               // so that L2 - L1 is the row count of the block X( K, L )

               L1 = (L - 1) * NB + 1
               L2 = MIN( L * NB, N ) + 1

               ztrsyl(TRANA, TRANB, ISGN, K2-K1, L2-L1, A( K1, K1 ), LDA, B( L1, L1 ), LDB, C( K1, L1 ), LDC, SCALOC, IINFO );
               INFO = MAX( INFO, IINFO )

               if ( SCALOC * SWORK( K, L ) .EQ. ZERO ) {
                  if ( SCALOC .EQ. ZERO ) {
                     // The magnitude of the largest entry of X(K1:K2-1, L1:L2-1)
                     // is larger than the product of BIGNUM**2 and cannot be
                     // represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1).
                     // Mark the computation as pointless.
                     BUF = ZERO
                  } else {
                     // Use second scaling factor to prevent flushing to zero.
                     BUF = BUF*2.D0**EXPONENT( SCALOC )
                  }
                  for (JJ = 1; JJ <= NBB; JJ++) {
                     for (LL = 1; LL <= NBA; LL++) {
                        // Bound by BIGNUM to not introduce Inf. The value
                        // is irrelevant; corresponding entries of the
                        // solution will be flushed in consistency scaling.
                        SWORK( LL, JJ ) = MIN( BIGNUM, SWORK( LL, JJ ) / 2.D0**EXPONENT( SCALOC ) )
                     END DO
                  END DO
               }
               SWORK( K, L ) = SCALOC * SWORK( K, L )
               XNRM = ZLANGE( 'I', K2-K1, L2-L1, C( K1, L1 ), LDC, WNRM )

               DO I = K + 1, NBA

                  // C( I, L ) := C( I, L ) - A( K, I )**H * C( K, L )

                  I1 = (I - 1) * NB + 1
                  I2 = MIN( I * NB, M ) + 1

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = ZLANGE( 'I', I2-I1, L2-L1, C( I1, L1 ), LDC, WNRM )
                  SCAMIN = MIN( SWORK( I, L ), SWORK( K, L ) )
                  CNRM = CNRM * ( SCAMIN / SWORK( I, L ) )
                  XNRM = XNRM * ( SCAMIN / SWORK( K, L ) )
                  ANRM = SWORK( I, AWRK + K )
                  SCALOC = DLARMM( ANRM, XNRM, CNRM )
                  if ( SCALOC * SCAMIN .EQ. ZERO ) {
                     // Use second scaling factor to prevent flushing to zero.
                     BUF = BUF*2.D0**EXPONENT( SCALOC )
                     for (JJ = 1; JJ <= NBB; JJ++) {
                        for (LL = 1; LL <= NBA; LL++) {
                        SWORK( LL, JJ ) = MIN( BIGNUM, SWORK( LL, JJ ) / 2.D0**EXPONENT( SCALOC ) )
                        END DO
                     END DO
                     SCAMIN = SCAMIN / 2.D0**EXPONENT( SCALOC )
                     SCALOC = SCALOC / 2.D0**EXPONENT( SCALOC )
                  }
                  CNRM = CNRM * SCALOC
                  XNRM = XNRM * SCALOC

                  // Simultaneously apply the robust update factor and the
                  // consistency scaling factor to to C( I, L ) and C( K, L).

                  SCAL = ( SCAMIN / SWORK( K, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     DO LL = L1, L2-1
                        zdscal(K2-K1, SCAL, C( K1, LL ), 1 );
                     END DO
                  ENDIF

                  SCAL = ( SCAMIN / SWORK( I, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     DO LL = L1, L2-1
                        zdscal(I2-I1, SCAL, C( I1, LL ), 1 );
                     END DO
                  ENDIF

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( I, L ) = SCAMIN * SCALOC

                  zgemm('C', 'N', I2-I1, L2-L1, K2-K1, -CONE, A( K1, I1 ), LDA, C( K1, L1 ), LDC, CONE, C( I1, L1 ), LDC );
               END DO

               DO J = L + 1, NBB

                  // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( L, J )

                  J1 = (J - 1) * NB + 1
                  J2 = MIN( J * NB, N ) + 1

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = ZLANGE( 'I', K2-K1, J2-J1, C( K1, J1 ), LDC, WNRM )
                  SCAMIN = MIN( SWORK( K, J ), SWORK( K, L ) )
                  CNRM = CNRM * ( SCAMIN / SWORK( K, J ) )
                  XNRM = XNRM * ( SCAMIN / SWORK( K, L ) )
                  BNRM = SWORK( L, BWRK + J )
                  SCALOC = DLARMM( BNRM, XNRM, CNRM )
                  if ( SCALOC * SCAMIN .EQ. ZERO ) {
                     // Use second scaling factor to prevent flushing to zero.
                     BUF = BUF*2.D0**EXPONENT( SCALOC )
                     for (JJ = 1; JJ <= NBB; JJ++) {
                        for (LL = 1; LL <= NBA; LL++) {
                        SWORK( LL, JJ ) = MIN( BIGNUM, SWORK( LL, JJ ) / 2.D0**EXPONENT( SCALOC ) )
                        END DO
                     END DO
                     SCAMIN = SCAMIN / 2.D0**EXPONENT( SCALOC )
                     SCALOC = SCALOC / 2.D0**EXPONENT( SCALOC )
                  }
                  CNRM = CNRM * SCALOC
                  XNRM = XNRM * SCALOC

                  // Simultaneously apply the robust update factor and the
                  // consistency scaling factor to to C( K, J ) and C( K, L).

                  SCAL = ( SCAMIN / SWORK( K, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                      DO LL = L1, L2-1
                         zdscal(K2-K1, SCAL, C( K1, LL ), 1 );
                      END DO
                  ENDIF

                  SCAL = ( SCAMIN / SWORK( K, J ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     DO JJ = J1, J2-1
                        zdscal(K2-K1, SCAL, C( K1, JJ ), 1 );
                     END DO
                  ENDIF

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( K, J ) = SCAMIN * SCALOC

                  zgemm('N', 'N', K2-K1, J2-J1, L2-L1, -CSGN, C( K1, L1 ), LDC, B( L1, J1 ), LDB, CONE, C( K1, J1 ), LDC );
               END DO
            END DO
         END DO
      } else if ( .NOT.NOTRNA .AND. .NOT.NOTRNB ) {

         // Solve    A**H *X + ISGN*X*B**H = scale*C.

         // The (K,L)th block of X is determined starting from
         // top-right corner column by column by

            // A(K,K)**H*X(K,L) + ISGN*X(K,L)*B(L,L)**H = C(K,L) - R(K,L)

         // Where
                      // K-1                          N
             // R(K,L) = SUM [A(I,K)**H*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**H].
                      // I=1                        J=L+1

         // Start loop over block rows (index = K) and block columns (index = L)

         for (K = 1; K <= NBA; K++) {

            // K1: row index of the first row in X( K, L )
            // K2: row index of the first row in X( K+1, L )
            // so the K2 - K1 is the column count of the block X( K, L )

            K1 = (K - 1) * NB + 1
            K2 = MIN( K * NB, M ) + 1
            DO L = NBB, 1, -1

               // L1: column index of the first column in X( K, L )
               // L2: column index of the first column in X( K, L + 1)
               // so that L2 - L1 is the row count of the block X( K, L )

               L1 = (L - 1) * NB + 1
               L2 = MIN( L * NB, N ) + 1

               ztrsyl(TRANA, TRANB, ISGN, K2-K1, L2-L1, A( K1, K1 ), LDA, B( L1, L1 ), LDB, C( K1, L1 ), LDC, SCALOC, IINFO );
               INFO = MAX( INFO, IINFO )

               if ( SCALOC * SWORK( K, L ) .EQ. ZERO ) {
                  if ( SCALOC .EQ. ZERO ) {
                     // The magnitude of the largest entry of X(K1:K2-1, L1:L2-1)
                     // is larger than the product of BIGNUM**2 and cannot be
                     // represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1).
                     // Mark the computation as pointless.
                     BUF = ZERO
                  } else {
                     // Use second scaling factor to prevent flushing to zero.
                     BUF = BUF*2.D0**EXPONENT( SCALOC )
                  }
                  for (JJ = 1; JJ <= NBB; JJ++) {
                     for (LL = 1; LL <= NBA; LL++) {
                        // Bound by BIGNUM to not introduce Inf. The value
                        // is irrelevant; corresponding entries of the
                        // solution will be flushed in consistency scaling.
                        SWORK( LL, JJ ) = MIN( BIGNUM, SWORK( LL, JJ ) / 2.D0**EXPONENT( SCALOC ) )
                     END DO
                  END DO
               }
               SWORK( K, L ) = SCALOC * SWORK( K, L )
               XNRM = ZLANGE( 'I', K2-K1, L2-L1, C( K1, L1 ), LDC, WNRM )

               DO I = K + 1, NBA

                  // C( I, L ) := C( I, L ) - A( K, I )**H * C( K, L )

                  I1 = (I - 1) * NB + 1
                  I2 = MIN( I * NB, M ) + 1

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = ZLANGE( 'I', I2-I1, L2-L1, C( I1, L1 ), LDC, WNRM )
                  SCAMIN = MIN( SWORK( I, L ), SWORK( K, L ) )
                  CNRM = CNRM * ( SCAMIN / SWORK( I, L ) )
                  XNRM = XNRM * ( SCAMIN / SWORK( K, L ) )
                  ANRM = SWORK( I, AWRK + K )
                  SCALOC = DLARMM( ANRM, XNRM, CNRM )
                  if ( SCALOC * SCAMIN .EQ. ZERO ) {
                     // Use second scaling factor to prevent flushing to zero.
                     BUF = BUF*2.D0**EXPONENT( SCALOC )
                     for (JJ = 1; JJ <= NBB; JJ++) {
                        for (LL = 1; LL <= NBA; LL++) {
                        SWORK( LL, JJ ) = MIN( BIGNUM, SWORK( LL, JJ ) / 2.D0**EXPONENT( SCALOC ) )
                        END DO
                     END DO
                     SCAMIN = SCAMIN / 2.D0**EXPONENT( SCALOC )
                     SCALOC = SCALOC / 2.D0**EXPONENT( SCALOC )
                  }
                  CNRM = CNRM * SCALOC
                  XNRM = XNRM * SCALOC

                  // Simultaneously apply the robust update factor and the
                  // consistency scaling factor to C( I, L ) and C( K, L).

                  SCAL = ( SCAMIN / SWORK( K, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     DO LL = L1, L2-1
                        zdscal(K2-K1, SCAL, C( K1, LL ), 1 );
                     END DO
                  ENDIF

                  SCAL = ( SCAMIN / SWORK( I, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     DO LL = L1, L2-1
                        zdscal(I2-I1, SCAL, C( I1, LL ), 1 );
                     END DO
                  ENDIF

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( I, L ) = SCAMIN * SCALOC

                  zgemm('C', 'N', I2-I1, L2-L1, K2-K1, -CONE, A( K1, I1 ), LDA, C( K1, L1 ), LDC, CONE, C( I1, L1 ), LDC );
               END DO

               DO J = 1, L - 1

                  // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( J, L )**H

                  J1 = (J - 1) * NB + 1
                  J2 = MIN( J * NB, N ) + 1

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = ZLANGE( 'I', K2-K1, J2-J1, C( K1, J1 ), LDC, WNRM )
                  SCAMIN = MIN( SWORK( K, J ), SWORK( K, L ) )
                  CNRM = CNRM * ( SCAMIN / SWORK( K, J ) )
                  XNRM = XNRM * ( SCAMIN / SWORK( K, L ) )
                  BNRM = SWORK( L, BWRK + J )
                  SCALOC = DLARMM( BNRM, XNRM, CNRM )
                  if ( SCALOC * SCAMIN .EQ. ZERO ) {
                     // Use second scaling factor to prevent flushing to zero.
                     BUF = BUF*2.D0**EXPONENT( SCALOC )
                     for (JJ = 1; JJ <= NBB; JJ++) {
                        for (LL = 1; LL <= NBA; LL++) {
                        SWORK( LL, JJ ) = MIN( BIGNUM, SWORK( LL, JJ ) / 2.D0**EXPONENT( SCALOC ) )
                        END DO
                     END DO
                     SCAMIN = SCAMIN / 2.D0**EXPONENT( SCALOC )
                     SCALOC = SCALOC / 2.D0**EXPONENT( SCALOC )
                  }
                  CNRM = CNRM * SCALOC
                  XNRM = XNRM * SCALOC

                  // Simultaneously apply the robust update factor and the
                  // consistency scaling factor to C( K, J ) and C( K, L).

                  SCAL = ( SCAMIN / SWORK( K, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     DO LL = L1, L2-1
                        zdscal(K2-K1, SCAL, C( K1, LL ), 1);
                     END DO
                  ENDIF

                  SCAL = ( SCAMIN / SWORK( K, J ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     DO JJ = J1, J2-1
                        zdscal(K2-K1, SCAL, C( K1, JJ ), 1 );
                     END DO
                  ENDIF

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( K, J ) = SCAMIN * SCALOC

                  zgemm('N', 'C', K2-K1, J2-J1, L2-L1, -CSGN, C( K1, L1 ), LDC, B( J1, L1 ), LDB, CONE, C( K1, J1 ), LDC );
               END DO
            END DO
         END DO
      } else if ( NOTRNA .AND. .NOT.NOTRNB ) {

         // Solve    A*X + ISGN*X*B**H = scale*C.

         // The (K,L)th block of X is determined starting from
         // bottom-right corner column by column by

             // A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)**H = C(K,L) - R(K,L)

         // Where
                       // M                          N
             // R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**H].
                     // I=K+1                      J=L+1

         // Start loop over block rows (index = K) and block columns (index = L)

         DO K = NBA, 1, -1

            // K1: row index of the first row in X( K, L )
            // K2: row index of the first row in X( K+1, L )
            // so the K2 - K1 is the column count of the block X( K, L )

            K1 = (K - 1) * NB + 1
            K2 = MIN( K * NB, M ) + 1
            DO L = NBB, 1, -1

               // L1: column index of the first column in X( K, L )
               // L2: column index of the first column in X( K, L + 1)
               // so that L2 - L1 is the row count of the block X( K, L )

               L1 = (L - 1) * NB + 1
               L2 = MIN( L * NB, N ) + 1

               ztrsyl(TRANA, TRANB, ISGN, K2-K1, L2-L1, A( K1, K1 ), LDA, B( L1, L1 ), LDB, C( K1, L1 ), LDC, SCALOC, IINFO );
               INFO = MAX( INFO, IINFO )

               if ( SCALOC * SWORK( K, L ) .EQ. ZERO ) {
                  if ( SCALOC .EQ. ZERO ) {
                     // The magnitude of the largest entry of X(K1:K2-1, L1:L2-1)
                     // is larger than the product of BIGNUM**2 and cannot be
                     // represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1).
                     // Mark the computation as pointless.
                     BUF = ZERO
                  } else {
                     // Use second scaling factor to prevent flushing to zero.
                     BUF = BUF*2.D0**EXPONENT( SCALOC )
                  }
                  for (JJ = 1; JJ <= NBB; JJ++) {
                     for (LL = 1; LL <= NBA; LL++) {
                        // Bound by BIGNUM to not introduce Inf. The value
                        // is irrelevant; corresponding entries of the
                        // solution will be flushed in consistency scaling.
                        SWORK( LL, JJ ) = MIN( BIGNUM, SWORK( LL, JJ ) / 2.D0**EXPONENT( SCALOC ) )
                     END DO
                  END DO
               }
               SWORK( K, L ) = SCALOC * SWORK( K, L )
               XNRM = ZLANGE( 'I', K2-K1, L2-L1, C( K1, L1 ), LDC, WNRM )

               DO I = 1, K - 1

                  // C( I, L ) := C( I, L ) - A( I, K ) * C( K, L )

                  I1 = (I - 1) * NB + 1
                  I2 = MIN( I * NB, M ) + 1

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = ZLANGE( 'I', I2-I1, L2-L1, C( I1, L1 ), LDC, WNRM )
                  SCAMIN = MIN( SWORK( I, L ), SWORK( K, L ) )
                  CNRM = CNRM * ( SCAMIN / SWORK( I, L ) )
                  XNRM = XNRM * ( SCAMIN / SWORK( K, L ) )
                  ANRM = SWORK( I, AWRK + K )
                  SCALOC = DLARMM( ANRM, XNRM, CNRM )
                  if ( SCALOC * SCAMIN .EQ. ZERO ) {
                     // Use second scaling factor to prevent flushing to zero.
                     BUF = BUF*2.D0**EXPONENT( SCALOC )
                     for (JJ = 1; JJ <= NBB; JJ++) {
                        for (LL = 1; LL <= NBA; LL++) {
                        SWORK( LL, JJ ) = MIN( BIGNUM, SWORK( LL, JJ ) / 2.D0**EXPONENT( SCALOC ) )
                        END DO
                     END DO
                     SCAMIN = SCAMIN / 2.D0**EXPONENT( SCALOC )
                     SCALOC = SCALOC / 2.D0**EXPONENT( SCALOC )
                  }
                  CNRM = CNRM * SCALOC
                  XNRM = XNRM * SCALOC

                  // Simultaneously apply the robust update factor and the
                  // consistency scaling factor to C( I, L ) and C( K, L).

                  SCAL = ( SCAMIN / SWORK( K, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     DO LL = L1, L2-1
                        zdscal(K2-K1, SCAL, C( K1, LL ), 1 );
                     END DO
                  ENDIF

                  SCAL = ( SCAMIN / SWORK( I, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     DO LL = L1, L2-1
                        zdscal(I2-I1, SCAL, C( I1, LL ), 1 );
                     END DO
                  ENDIF

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( I, L ) = SCAMIN * SCALOC

                  zgemm('N', 'N', I2-I1, L2-L1, K2-K1, -CONE, A( I1, K1 ), LDA, C( K1, L1 ), LDC, CONE, C( I1, L1 ), LDC );

               END DO

               DO J = 1, L - 1

                  // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( J, L )**H

                  J1 = (J - 1) * NB + 1
                  J2 = MIN( J * NB, N ) + 1

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = ZLANGE( 'I', K2-K1, J2-J1, C( K1, J1 ), LDC, WNRM )
                  SCAMIN = MIN( SWORK( K, J ), SWORK( K, L ) )
                  CNRM = CNRM * ( SCAMIN / SWORK( K, J ) )
                  XNRM = XNRM * ( SCAMIN / SWORK( K, L ) )
                  BNRM = SWORK( L, BWRK + J )
                  SCALOC = DLARMM( BNRM, XNRM, CNRM )
                  if ( SCALOC * SCAMIN .EQ. ZERO ) {
                     // Use second scaling factor to prevent flushing to zero.
                     BUF = BUF*2.D0**EXPONENT( SCALOC )
                     for (JJ = 1; JJ <= NBB; JJ++) {
                        for (LL = 1; LL <= NBA; LL++) {
                        SWORK( LL, JJ ) = MIN( BIGNUM, SWORK( LL, JJ ) / 2.D0**EXPONENT( SCALOC ) )
                        END DO
                     END DO
                     SCAMIN = SCAMIN / 2.D0**EXPONENT( SCALOC )
                     SCALOC = SCALOC / 2.D0**EXPONENT( SCALOC )
                  }
                  CNRM = CNRM * SCALOC
                  XNRM = XNRM * SCALOC

                  // Simultaneously apply the robust update factor and the
                  // consistency scaling factor to C( K, J ) and C( K, L).

                  SCAL = ( SCAMIN / SWORK( K, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     DO JJ = L1, L2-1
                        zdscal(K2-K1, SCAL, C( K1, JJ ), 1 );
                     END DO
                  ENDIF

                  SCAL = ( SCAMIN / SWORK( K, J ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     DO JJ = J1, J2-1
                        zdscal(K2-K1, SCAL, C( K1, JJ ), 1 );
                     END DO
                  ENDIF

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( K, J ) = SCAMIN * SCALOC

                  zgemm('N', 'C', K2-K1, J2-J1, L2-L1, -CSGN, C( K1, L1 ), LDC, B( J1, L1 ), LDB, CONE, C( K1, J1 ), LDC );
               END DO
            END DO
         END DO

      }

      // Reduce local scaling factors

      SCALE = SWORK( 1, 1 )
      for (K = 1; K <= NBA; K++) {
         for (L = 1; L <= NBB; L++) {
            SCALE = MIN( SCALE, SWORK( K, L ) )
         END DO
      END DO
      if ( SCALE .EQ. ZERO ) {

         // The magnitude of the largest entry of the solution is larger
         // than the product of BIGNUM**2 and cannot be represented in the
         // form (1/SCALE)*X if SCALE is double          . Set SCALE to;
         // zero and give up.

         SWORK(1,1) = MAX( NBA, NBB )
         SWORK(2,1) = 2 * NBB + NBA
         RETURN
      }

      // Realize consistent scaling

      for (K = 1; K <= NBA; K++) {
         K1 = (K - 1) * NB + 1
         K2 = MIN( K * NB, M ) + 1
         for (L = 1; L <= NBB; L++) {
            L1 = (L - 1) * NB + 1
            L2 = MIN( L * NB, N ) + 1
            SCAL = SCALE / SWORK( K, L )
            if ( SCAL .NE. ONE ) {
               DO LL = L1, L2-1
                  zdscal(K2-K1, SCAL, C( K1, LL ), 1 );
               END DO
            ENDIF
         END DO
      END DO

      if ( BUF .NE. ONE .AND. BUF.GT.ZERO ) {

         // Decrease SCALE as much as possible.

         SCALOC = MIN( SCALE / SMLNUM, ONE / BUF )
         BUF = BUF * SCALOC
         SCALE = SCALE / SCALOC
      }

      if ( BUF.NE.ONE .AND. BUF.GT.ZERO ) {

         // In case of overly aggressive scaling during the computation,
         // flushing of the global scale factor may be prevented by
         // undoing some of the scaling. This step is to ensure that
         // this routine flushes only scale factors that TRSYL also
         // flushes and be usable as a drop-in replacement.

         // How much can the normwise largest entry be upscaled?

         SCAL = MAX( ABS( DBLE( C( 1, 1 ) ) ), ABS( DIMAG( C ( 1, 1 ) ) ) )
         for (K = 1; K <= M; K++) {
            for (L = 1; L <= N; L++) {
               SCAL = MAX( SCAL, ABS( DBLE ( C( K, L ) ) ), ABS( DIMAG ( C( K, L ) ) ) )
            END DO
         END DO

         // Increase BUF as close to 1 as possible and apply scaling.

         SCALOC = MIN( BIGNUM / SCAL, ONE / BUF )
         BUF = BUF * SCALOC
         zlascl('G', -1, -1, ONE, SCALOC, M, N, C, LDC, IINFO );
      }

      // Combine with buffer scaling factor. SCALE will be flushed if
      // BUF is less than one here.

      SCALE = SCALE * BUF

      // Restore workspace dimensions

      SWORK(1,1) = MAX( NBA, NBB )
      SWORK(2,1) = 2 * NBB + NBA

      RETURN

      // End of ZTRSYL3

      }
