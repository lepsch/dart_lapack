      SUBROUTINE DTRSYL3( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, LDC, SCALE, IWORK, LIWORK, SWORK, LDSWORK, INFO )
      IMPLICIT NONE

      // .. Scalar Arguments ..
      String             TRANA, TRANB;
      int                INFO, ISGN, LDA, LDB, LDC, M, N, LIWORK, LDSWORK;
      double             SCALE;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             A( LDA, * ), B( LDB, * ), C( LDC, * ), SWORK( LDSWORK, * );
      // ..
      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRNA, NOTRNB, LQUERY, SKIP;
      int                AWRK, BWRK, I, I1, I2, IINFO, J, J1, J2, JJ, K, K1, K2, L, L1, L2, LL, NBA, NB, NBB, PC;
      double             ANRM, BIGNUM, BNRM, CNRM, SCAL, SCALOC, SCAMIN, SGN, XNRM, BUF, SMLNUM;
      // ..
      // .. Local Arrays ..
      double             WNRM( MAX( M, N ) );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLANGE, DLAMCH, DLARMM;
      // EXTERNAL DLANGE, DLAMCH, DLARMM, ILAENV, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLASCL, DSCAL, DTRSYL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, EXPONENT, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Decode and Test input parameters

      NOTRNA = LSAME( TRANA, 'N' )
      NOTRNB = LSAME( TRANB, 'N' )

      // Use the same block size for all matrices.

      NB = MAX(8, ILAENV( 1, 'DTRSYL', '', M, N, -1, -1) )

      // Compute number of blocks in A and B

      NBA = MAX( 1, (M + NB - 1) / NB )
      NBB = MAX( 1, (N + NB - 1) / NB )

      // Compute workspace

      INFO = 0
      LQUERY = ( LIWORK.EQ.-1 .OR. LDSWORK.EQ.-1 )
      IWORK( 1 ) = NBA + NBB + 2
      if ( LQUERY ) {
         LDSWORK = 2
         SWORK( 1, 1 ) = MAX( NBA, NBB )
         SWORK( 2, 1 ) = 2 * NBB + NBA
      }

      // Test the input arguments

      if ( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'T' ) .AND. .NOT. LSAME( TRANA, 'C' ) ) {
         INFO = -1
      } else if ( .NOT.NOTRNB .AND. .NOT.LSAME( TRANB, 'T' ) .AND. .NOT. LSAME( TRANB, 'C' ) ) {
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
         xerbla('DTRSYL3', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      SCALE = ONE
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      // Use unblocked code for small problems or if insufficient
      // workspaces are provided

      IF( MIN( NBA, NBB ).EQ.1 .OR. LDSWORK.LT.MAX( NBA, NBB ) .OR. LIWORK.LT.IWORK(1) ) THEN         CALL DTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, LDC, SCALE, INFO )
        RETURN
      }

      // Set constants to control overflow

      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM

       // Partition A such that 2-by-2 blocks on the diagonal are not split

       SKIP = .FALSE.
       for (I = 1; I <= NBA; I++) {
          IWORK( I ) = ( I - 1 ) * NB + 1
       END DO
       IWORK( NBA + 1 ) = M + 1
       for (K = 1; K <= NBA; K++) {
          L1 = IWORK( K )
          L2 = IWORK( K + 1 ) - 1
          for (L = L1; L <= L2; L++) {
             if ( SKIP ) {
                SKIP = .FALSE.
                CYCLE
             }
             if ( L.GE.M ) {
                // A( M, M ) is a 1-by-1 block
                CYCLE
             }
             if ( A( L, L+1 ).NE.ZERO .AND. A( L+1, L ).NE.ZERO ) {
                // Check if 2-by-2 block is split
                if ( L + 1 .EQ. IWORK( K + 1 ) ) {
                   IWORK( K + 1 ) = IWORK( K + 1 ) + 1
                   CYCLE
                }
                SKIP = .TRUE.
             }
          END DO
       END DO
       IWORK( NBA + 1 ) = M + 1
       if ( IWORK( NBA ).GE.IWORK( NBA + 1 ) ) {
          IWORK( NBA ) = IWORK( NBA + 1 )
          NBA = NBA - 1
       }

       // Partition B such that 2-by-2 blocks on the diagonal are not split

       PC = NBA + 1
       SKIP = .FALSE.
       for (I = 1; I <= NBB; I++) {
          IWORK( PC + I ) = ( I - 1 ) * NB + 1
       END DO
       IWORK( PC + NBB + 1 ) = N + 1
       for (K = 1; K <= NBB; K++) {
          L1 = IWORK( PC + K )
          L2 = IWORK( PC + K + 1 ) - 1
          for (L = L1; L <= L2; L++) {
             if ( SKIP ) {
                SKIP = .FALSE.
                CYCLE
             }
             if ( L.GE.N ) {
                // B( N, N ) is a 1-by-1 block
                CYCLE
             }
             if ( B( L, L+1 ).NE.ZERO .AND. B( L+1, L ).NE.ZERO ) {
                // Check if 2-by-2 block is split
                if ( L + 1 .EQ. IWORK( PC + K + 1 ) ) {
                   IWORK( PC + K + 1 ) = IWORK( PC + K + 1 ) + 1
                   CYCLE
                }
                SKIP = .TRUE.
             }
          END DO
       END DO
       IWORK( PC + NBB + 1 ) = N + 1
       if ( IWORK( PC + NBB ).GE.IWORK( PC + NBB + 1 ) ) {
          IWORK( PC + NBB ) = IWORK( PC + NBB + 1 )
          NBB = NBB - 1
       }

      // Set local scaling factors - must never attain zero.

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
         K1 = IWORK( K )
         K2 = IWORK( K + 1 )
         for (L = K; L <= NBA; L++) {
            L1 = IWORK( L )
            L2 = IWORK( L + 1 )
            if ( NOTRNA ) {
               SWORK( K, AWRK + L ) = DLANGE( 'I', K2-K1, L2-L1, A( K1, L1 ), LDA, WNRM )
            } else {
               SWORK( L, AWRK + K ) = DLANGE( '1', K2-K1, L2-L1, A( K1, L1 ), LDA, WNRM )
            }
         END DO
      END DO
      BWRK = NBB + NBA
      for (K = 1; K <= NBB; K++) {
         K1 = IWORK( PC + K )
         K2 = IWORK( PC + K + 1 )
         for (L = K; L <= NBB; L++) {
            L1 = IWORK( PC + L )
            L2 = IWORK( PC + L + 1 )
            if ( NOTRNB ) {
               SWORK( K, BWRK + L ) = DLANGE( 'I', K2-K1, L2-L1, B( K1, L1 ), LDB, WNRM )
            } else {
               SWORK( L, BWRK + K ) = DLANGE( '1', K2-K1, L2-L1, B( K1, L1 ), LDB, WNRM )
            }
         END DO
      END DO

      SGN = DBLE( ISGN )

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

            K1 = IWORK( K )
            K2 = IWORK( K + 1 )
            for (L = 1; L <= NBB; L++) {

               // L1: column index of the first column in X( K, L )
               // L2: column index of the first column in X( K, L + 1)
               // so that L2 - L1 is the row count of the block X( K, L )

               L1 = IWORK( PC + L )
               L2 = IWORK( PC + L + 1 )

               dtrsyl(TRANA, TRANB, ISGN, K2-K1, L2-L1, A( K1, K1 ), LDA, B( L1, L1 ), LDB, C( K1, L1 ), LDC, SCALOC, IINFO );
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
               XNRM = DLANGE( 'I', K2-K1, L2-L1, C( K1, L1 ), LDC, WNRM )

               DO I = K - 1, 1, -1

                  // C( I, L ) := C( I, L ) - A( I, K ) * C( K, L )

                  I1 = IWORK( I )
                  I2 = IWORK( I + 1 )

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = DLANGE( 'I', I2-I1, L2-L1, C( I1, L1 ), LDC, WNRM )
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
                  if (SCAL .NE. ONE) {
                      for (JJ = L1; JJ <= L2-1; JJ++) {
                         dscal(K2-K1, SCAL, C( K1, JJ ), 1);
                      END DO
                  }

                  SCAL = ( SCAMIN / SWORK( I, L ) ) * SCALOC
                  if (SCAL .NE. ONE) {
                      for (LL = L1; LL <= L2-1; LL++) {
                         dscal(I2-I1, SCAL, C( I1, LL ), 1);
                      END DO
                  }

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( I, L ) = SCAMIN * SCALOC

                  dgemm('N', 'N', I2-I1, L2-L1, K2-K1, -ONE, A( I1, K1 ), LDA, C( K1, L1 ), LDC, ONE, C( I1, L1 ), LDC );

               END DO

               for (J = L + 1; J <= NBB; J++) {

                  // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( L, J )

                  J1 = IWORK( PC + J )
                  J2 = IWORK( PC + J + 1 )

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = DLANGE( 'I', K2-K1, J2-J1, C( K1, J1 ), LDC, WNRM )
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
                     for (LL = L1; LL <= L2-1; LL++) {
                        dscal(K2-K1, SCAL, C( K1, LL ), 1 );
                     END DO
                  }

                  SCAL = ( SCAMIN / SWORK( K, J ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                      for (JJ = J1; JJ <= J2-1; JJ++) {
                         dscal(K2-K1, SCAL, C( K1, JJ ), 1 );
                      END DO
                  }

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( K, J ) = SCAMIN * SCALOC

                  dgemm('N', 'N', K2-K1, J2-J1, L2-L1, -SGN, C( K1, L1 ), LDC, B( L1, J1 ), LDB, ONE, C( K1, J1 ), LDC );
               END DO
            END DO
         END DO
      } else if ( .NOT.NOTRNA .AND. NOTRNB ) {

         // Solve    A**T*X + ISGN*X*B = scale*C.

         // The (K,L)th block of X is determined starting from
         // upper-left corner column by column by

           // A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

         // Where
                    // K-1                        L-1
           // R(K,L) = SUM [A(I,K)**T*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
                    // I=1                        J=1

         // Start loop over block rows (index = K) and block columns (index = L)

         for (K = 1; K <= NBA; K++) {

            // K1: row index of the first row in X( K, L )
            // K2: row index of the first row in X( K+1, L )
            // so the K2 - K1 is the column count of the block X( K, L )

            K1 = IWORK( K )
            K2 = IWORK( K + 1 )
            for (L = 1; L <= NBB; L++) {

               // L1: column index of the first column in X( K, L )
               // L2: column index of the first column in X( K, L + 1)
               // so that L2 - L1 is the row count of the block X( K, L )

               L1 = IWORK( PC + L )
               L2 = IWORK( PC + L + 1 )

               dtrsyl(TRANA, TRANB, ISGN, K2-K1, L2-L1, A( K1, K1 ), LDA, B( L1, L1 ), LDB, C( K1, L1 ), LDC, SCALOC, IINFO );
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
               XNRM = DLANGE( 'I', K2-K1, L2-L1, C( K1, L1 ), LDC, WNRM )

               for (I = K + 1; I <= NBA; I++) {

                  // C( I, L ) := C( I, L ) - A( K, I )**T * C( K, L )

                  I1 = IWORK( I )
                  I2 = IWORK( I + 1 )

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = DLANGE( 'I', I2-I1, L2-L1, C( I1, L1 ), LDC, WNRM )
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
                  // consistency scaling factor to to C( I, L ) and C( K, L ).

                  SCAL = ( SCAMIN / SWORK( K, L ) ) * SCALOC
                  if (SCAL .NE. ONE) {
                     for (LL = L1; LL <= L2-1; LL++) {
                        dscal(K2-K1, SCAL, C( K1, LL ), 1 );
                     END DO
                  }

                  SCAL = ( SCAMIN / SWORK( I, L ) ) * SCALOC
                  if (SCAL .NE. ONE) {
                     for (LL = L1; LL <= L2-1; LL++) {
                        dscal(I2-I1, SCAL, C( I1, LL ), 1 );
                     END DO
                  }

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( I, L ) = SCAMIN * SCALOC

                  dgemm('T', 'N', I2-I1, L2-L1, K2-K1, -ONE, A( K1, I1 ), LDA, C( K1, L1 ), LDC, ONE, C( I1, L1 ), LDC );
               END DO

               for (J = L + 1; J <= NBB; J++) {

                  // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( L, J )

                  J1 = IWORK( PC + J )
                  J2 = IWORK( PC + J + 1 )

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = DLANGE( 'I', K2-K1, J2-J1, C( K1, J1 ), LDC, WNRM )
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
                  // consistency scaling factor to to C( K, J ) and C( K, L ).

                  SCAL = ( SCAMIN / SWORK( K, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                      for (LL = L1; LL <= L2-1; LL++) {
                         dscal(K2-K1, SCAL, C( K1, LL ), 1 );
                      END DO
                  }

                  SCAL = ( SCAMIN / SWORK( K, J ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     for (JJ = J1; JJ <= J2-1; JJ++) {
                        dscal(K2-K1, SCAL, C( K1, JJ ), 1 );
                     END DO
                  }

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( K, J ) = SCAMIN * SCALOC

                  dgemm('N', 'N', K2-K1, J2-J1, L2-L1, -SGN, C( K1, L1 ), LDC, B( L1, J1 ), LDB, ONE, C( K1, J1 ), LDC );
               END DO
            END DO
         END DO
      } else if ( .NOT.NOTRNA .AND. .NOT.NOTRNB ) {

         // Solve    A**T*X + ISGN*X*B**T = scale*C.

         // The (K,L)th block of X is determined starting from
         // top-right corner column by column by

            // A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)

         // Where
                      // K-1                          N
             // R(K,L) = SUM [A(I,K)**T*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
                      // I=1                        J=L+1

         // Start loop over block rows (index = K) and block columns (index = L)

         for (K = 1; K <= NBA; K++) {

            // K1: row index of the first row in X( K, L )
            // K2: row index of the first row in X( K+1, L )
            // so the K2 - K1 is the column count of the block X( K, L )

            K1 = IWORK( K )
            K2 = IWORK( K + 1 )
            DO L = NBB, 1, -1

               // L1: column index of the first column in X( K, L )
               // L2: column index of the first column in X( K, L + 1)
               // so that L2 - L1 is the row count of the block X( K, L )

               L1 = IWORK( PC + L )
               L2 = IWORK( PC + L + 1 )

               dtrsyl(TRANA, TRANB, ISGN, K2-K1, L2-L1, A( K1, K1 ), LDA, B( L1, L1 ), LDB, C( K1, L1 ), LDC, SCALOC, IINFO );
               INFO = MAX( INFO, IINFO )

               SWORK( K, L ) = SCALOC * SWORK( K, L )
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
               XNRM = DLANGE( 'I', K2-K1, L2-L1, C( K1, L1 ), LDC, WNRM )

               for (I = K + 1; I <= NBA; I++) {

                  // C( I, L ) := C( I, L ) - A( K, I )**T * C( K, L )

                  I1 = IWORK( I )
                  I2 = IWORK( I + 1 )

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = DLANGE( 'I', I2-I1, L2-L1, C( I1, L1 ), LDC, WNRM )
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
                  if (SCAL .NE. ONE) {
                     for (LL = L1; LL <= L2-1; LL++) {
                        dscal(K2-K1, SCAL, C( K1, LL ), 1 );
                     END DO
                  }

                  SCAL = ( SCAMIN / SWORK( I, L ) ) * SCALOC
                  if (SCAL .NE. ONE) {
                     for (LL = L1; LL <= L2-1; LL++) {
                        dscal(I2-I1, SCAL, C( I1, LL ), 1 );
                     END DO
                  }

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( I, L ) = SCAMIN * SCALOC

                  dgemm('T', 'N', I2-I1, L2-L1, K2-K1, -ONE, A( K1, I1 ), LDA, C( K1, L1 ), LDC, ONE, C( I1, L1 ), LDC );
               END DO

               for (J = 1; J <= L - 1; J++) {

                  // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( J, L )**T

                  J1 = IWORK( PC + J )
                  J2 = IWORK( PC + J + 1 )

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = DLANGE( 'I', K2-K1, J2-J1, C( K1, J1 ), LDC, WNRM )
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
                  // consistency scaling factor to C( K, J ) and C( K, L ).

                  SCAL = ( SCAMIN / SWORK( K, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     for (LL = L1; LL <= L2-1; LL++) {
                        dscal(K2-K1, SCAL, C( K1, LL ), 1);
                     END DO
                  }

                  SCAL = ( SCAMIN / SWORK( K, J ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     for (JJ = J1; JJ <= J2-1; JJ++) {
                        dscal(K2-K1, SCAL, C( K1, JJ ), 1 );
                     END DO
                  }

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( K, J ) = SCAMIN * SCALOC

                  dgemm('N', 'T', K2-K1, J2-J1, L2-L1, -SGN, C( K1, L1 ), LDC, B( J1, L1 ), LDB, ONE, C( K1, J1 ), LDC );
               END DO
            END DO
         END DO
      } else if ( NOTRNA .AND. .NOT.NOTRNB ) {

         // Solve    A*X + ISGN*X*B**T = scale*C.

         // The (K,L)th block of X is determined starting from
         // bottom-right corner column by column by

             // A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)

         // Where
                       // M                          N
             // R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
                     // I=K+1                      J=L+1

         // Start loop over block rows (index = K) and block columns (index = L)

         DO K = NBA, 1, -1

            // K1: row index of the first row in X( K, L )
            // K2: row index of the first row in X( K+1, L )
            // so the K2 - K1 is the column count of the block X( K, L )

            K1 = IWORK( K )
            K2 = IWORK( K + 1 )
            DO L = NBB, 1, -1

               // L1: column index of the first column in X( K, L )
               // L2: column index of the first column in X( K, L + 1)
               // so that L2 - L1 is the row count of the block X( K, L )

               L1 = IWORK( PC + L )
               L2 = IWORK( PC + L + 1 )

               dtrsyl(TRANA, TRANB, ISGN, K2-K1, L2-L1, A( K1, K1 ), LDA, B( L1, L1 ), LDB, C( K1, L1 ), LDC, SCALOC, IINFO );
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
               XNRM = DLANGE( 'I', K2-K1, L2-L1, C( K1, L1 ), LDC, WNRM )

               for (I = 1; I <= K - 1; I++) {

                  // C( I, L ) := C( I, L ) - A( I, K ) * C( K, L )

                  I1 = IWORK( I )
                  I2 = IWORK( I + 1 )

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = DLANGE( 'I', I2-I1, L2-L1, C( I1, L1 ), LDC, WNRM )
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
                  if (SCAL .NE. ONE) {
                     for (LL = L1; LL <= L2-1; LL++) {
                        dscal(K2-K1, SCAL, C( K1, LL ), 1 );
                     END DO
                  }

                  SCAL = ( SCAMIN / SWORK( I, L ) ) * SCALOC
                  if (SCAL .NE. ONE) {
                     for (LL = L1; LL <= L2-1; LL++) {
                        dscal(I2-I1, SCAL, C( I1, LL ), 1 );
                     END DO
                  }

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( I, L ) = SCAMIN * SCALOC

                  dgemm('N', 'N', I2-I1, L2-L1, K2-K1, -ONE, A( I1, K1 ), LDA, C( K1, L1 ), LDC, ONE, C( I1, L1 ), LDC );

               END DO

               for (J = 1; J <= L - 1; J++) {

                  // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( J, L )**T

                  J1 = IWORK( PC + J )
                  J2 = IWORK( PC + J + 1 )

                  // Compute scaling factor to survive the linear update
                  // simulating consistent scaling.

                  CNRM = DLANGE( 'I', K2-K1, J2-J1, C( K1, J1 ), LDC, WNRM )
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
                  // consistency scaling factor to C( K, J ) and C( K, L ).

                  SCAL = ( SCAMIN / SWORK( K, L ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     for (JJ = L1; JJ <= L2-1; JJ++) {
                        dscal(K2-K1, SCAL, C( K1, JJ ), 1 );
                     END DO
                  }

                  SCAL = ( SCAMIN / SWORK( K, J ) ) * SCALOC
                  if ( SCAL .NE. ONE ) {
                     for (JJ = J1; JJ <= J2-1; JJ++) {
                        dscal(K2-K1, SCAL, C( K1, JJ ), 1 );
                     END DO
                  }

                  // Record current scaling factor

                  SWORK( K, L ) = SCAMIN * SCALOC
                  SWORK( K, J ) = SCAMIN * SCALOC

                  dgemm('N', 'T', K2-K1, J2-J1, L2-L1, -SGN, C( K1, L1 ), LDC, B( J1, L1 ), LDB, ONE, C( K1, J1 ), LDC );
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

         IWORK(1) = NBA + NBB + 2
         SWORK(1,1) = MAX( NBA, NBB )
         SWORK(2,1) = 2 * NBB + NBA
         RETURN
      }

      // Realize consistent scaling

      for (K = 1; K <= NBA; K++) {
         K1 = IWORK( K )
         K2 = IWORK( K + 1 )
         for (L = 1; L <= NBB; L++) {
            L1 = IWORK( PC + L )
            L2 = IWORK( PC + L + 1 )
            SCAL = SCALE / SWORK( K, L )
            if ( SCAL .NE. ONE ) {
               for (LL = L1; LL <= L2-1; LL++) {
                  dscal(K2-K1, SCAL, C( K1, LL ), 1 );
               END DO
            }
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

         SCAL = C( 1, 1 )
         for (K = 1; K <= M; K++) {
            for (L = 1; L <= N; L++) {
               SCAL = MAX( SCAL, ABS( C( K, L ) ) )
            END DO
         END DO

         // Increase BUF as close to 1 as possible and apply scaling.

         SCALOC = MIN( BIGNUM / SCAL, ONE / BUF )
         BUF = BUF * SCALOC
         dlascl('G', -1, -1, ONE, SCALOC, M, N, C, LDC, IWORK(1) );
      }

      // Combine with buffer scaling factor. SCALE will be flushed if
      // BUF is less than one here.

      SCALE = SCALE * BUF

      // Restore workspace dimensions

      IWORK(1) = NBA + NBB + 2
      SWORK(1,1) = MAX( NBA, NBB )
      SWORK(2,1) = 2 * NBB + NBA

      RETURN

      // End of DTRSYL3

      }
