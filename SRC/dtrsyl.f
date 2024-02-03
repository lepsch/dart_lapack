      SUBROUTINE DTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, LDC, SCALE, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANA, TRANB;
      int                INFO, ISGN, LDA, LDB, LDC, M, N;
      double             SCALE;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), C( LDC, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRNA, NOTRNB;
      int                IERR, J, K, K1, K2, KNEXT, L, L1, L2, LNEXT;
      double             A11, BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN, SMLNUM, SUML, SUMR, XNORM;
      // ..
      // .. Local Arrays ..
      double             DUM( 1 ), VEC( 2, 2 ), X( 2, 2 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DDOT, DLAMCH, DLANGE;
      // EXTERNAL LSAME, DDOT, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLALN2, DLASY2, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Decode and Test input parameters

      NOTRNA = LSAME( TRANA, 'N' )
      NOTRNB = LSAME( TRANB, 'N' )

      INFO = 0
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
         xerbla('DTRSYL', -INFO );
         RETURN
      }

      // Quick return if possible

      SCALE = ONE
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      // Set constants to control overflow

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SMLNUM*DBLE( M*N ) / EPS
      BIGNUM = ONE / SMLNUM

      SMIN = MAX( SMLNUM, EPS*DLANGE( 'M', M, M, A, LDA, DUM ), EPS*DLANGE( 'M', N, N, B, LDB, DUM ) )

      SGN = ISGN

      if ( NOTRNA .AND. NOTRNB ) {

         // Solve    A*X + ISGN*X*B = scale*C.

         // The (K,L)th block of X is determined starting from
         // bottom-left corner column by column by

          // A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

         // Where
                   // M                         L-1
         // R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
                 // I=K+1                       J=1

         // Start column loop (index = L)
         // L1 (L2) : column index of the first (first) row of X(K,L).

         LNEXT = 1
         DO 60 L = 1, N
            IF( L.LT.LNEXT ) GO TO 60
            if ( L.EQ.N ) {
               L1 = L
               L2 = L
            } else {
               if ( B( L+1, L ).NE.ZERO ) {
                  L1 = L
                  L2 = L + 1
                  LNEXT = L + 2
               } else {
                  L1 = L
                  L2 = L
                  LNEXT = L + 1
               }
            }

            // Start row loop (index = K)
            // K1 (K2): row index of the first (last) row of X(K,L).

            KNEXT = M
            DO 50 K = M, 1, -1
               IF( K.GT.KNEXT ) GO TO 50
               if ( K.EQ.1 ) {
                  K1 = K
                  K2 = K
               } else {
                  if ( A( K, K-1 ).NE.ZERO ) {
                     K1 = K - 1
                     K2 = K
                     KNEXT = K - 2
                  } else {
                     K1 = K
                     K2 = K
                     KNEXT = K - 1
                  }
               }

               if ( L1.EQ.L2 .AND. K1.EQ.K2 ) {
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE

                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  if ( DA11.LE.SMIN ) {
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  }
                  DB = ABS( VEC( 1, 1 ) )
                  if ( DA11.LT.ONE .AND. DB.GT.ONE ) {
                     IF( DB.GT.BIGNUM*DA11 ) SCALOC = ONE / DB
                  }
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11

                  if ( SCALOC.NE.ONE ) {
                     DO 10 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
   10                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )

               } else if ( L1.EQ.L2 .AND. K1.NE.K2 ) {

                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )

                  dlaln2(.FALSE., 2, 1, SMIN, ONE, A( K1, K1 ), LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  IF( IERR.NE.0 ) INFO = 1

                  if ( SCALOC.NE.ONE ) {
                     DO 20 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
   20                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )

               } else if ( L1.NE.L2 .AND. K1.EQ.K2 ) {

                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )

                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, C( MIN( K1+1, M ), L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )

                  dlaln2(.TRUE., 2, 1, SMIN, ONE, B( L1, L1 ), LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  IF( IERR.NE.0 ) INFO = 1

                  if ( SCALOC.NE.ONE ) {
                     DO 30 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
   30                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )

               } else if ( L1.NE.L2 .AND. K1.NE.K2 ) {

                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )

                  dlasy2(.FALSE., .FALSE., ISGN, 2, 2, A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, 2, XNORM, IERR );
                  IF( IERR.NE.0 ) INFO = 1

                  if ( SCALOC.NE.ONE ) {
                     DO 40 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
   40                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               }

   50       CONTINUE

   60    CONTINUE

      } else if ( .NOT.NOTRNA .AND. NOTRNB ) {

         // Solve    A**T *X + ISGN*X*B = scale*C.

         // The (K,L)th block of X is determined starting from
         // upper-left corner column by column by

           // A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

         // Where
                    // K-1        T                    L-1
           // R(K,L) = SUM [A(I,K)**T*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
                    // I=1                          J=1

         // Start column loop (index = L)
         // L1 (L2): column index of the first (last) row of X(K,L)

         LNEXT = 1
         DO 120 L = 1, N
            IF( L.LT.LNEXT ) GO TO 120
            if ( L.EQ.N ) {
               L1 = L
               L2 = L
            } else {
               if ( B( L+1, L ).NE.ZERO ) {
                  L1 = L
                  L2 = L + 1
                  LNEXT = L + 2
               } else {
                  L1 = L
                  L2 = L
                  LNEXT = L + 1
               }
            }

            // Start row loop (index = K)
            // K1 (K2): row index of the first (last) row of X(K,L)

            KNEXT = 1
            DO 110 K = 1, M
               IF( K.LT.KNEXT ) GO TO 110
               if ( K.EQ.M ) {
                  K1 = K
                  K2 = K
               } else {
                  if ( A( K+1, K ).NE.ZERO ) {
                     K1 = K
                     K2 = K + 1
                     KNEXT = K + 2
                  } else {
                     K1 = K
                     K2 = K
                     KNEXT = K + 1
                  }
               }

               if ( L1.EQ.L2 .AND. K1.EQ.K2 ) {
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE

                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  if ( DA11.LE.SMIN ) {
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  }
                  DB = ABS( VEC( 1, 1 ) )
                  if ( DA11.LT.ONE .AND. DB.GT.ONE ) {
                     IF( DB.GT.BIGNUM*DA11 ) SCALOC = ONE / DB
                  }
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11

                  if ( SCALOC.NE.ONE ) {
                     DO 70 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
   70                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )

               } else if ( L1.EQ.L2 .AND. K1.NE.K2 ) {

                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )

                  dlaln2(.TRUE., 2, 1, SMIN, ONE, A( K1, K1 ), LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  IF( IERR.NE.0 ) INFO = 1

                  if ( SCALOC.NE.ONE ) {
                     DO 80 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
   80                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )

               } else if ( L1.NE.L2 .AND. K1.EQ.K2 ) {

                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )

                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )

                  dlaln2(.TRUE., 2, 1, SMIN, ONE, B( L1, L1 ), LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  IF( IERR.NE.0 ) INFO = 1

                  if ( SCALOC.NE.ONE ) {
                     DO 90 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
   90                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )

               } else if ( L1.NE.L2 .AND. K1.NE.K2 ) {

                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )

                  dlasy2(.TRUE., .FALSE., ISGN, 2, 2, A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, 2, XNORM, IERR );
                  IF( IERR.NE.0 ) INFO = 1

                  if ( SCALOC.NE.ONE ) {
                     DO 100 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
  100                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               }

  110       CONTINUE
  120    CONTINUE

      } else if ( .NOT.NOTRNA .AND. .NOT.NOTRNB ) {

         // Solve    A**T*X + ISGN*X*B**T = scale*C.

         // The (K,L)th block of X is determined starting from
         // top-right corner column by column by

            // A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)

         // Where
                      // K-1                            N
             // R(K,L) = SUM [A(I,K)**T*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
                      // I=1                          J=L+1

         // Start column loop (index = L)
         // L1 (L2): column index of the first (last) row of X(K,L)

         LNEXT = N
         DO 180 L = N, 1, -1
            IF( L.GT.LNEXT ) GO TO 180
            if ( L.EQ.1 ) {
               L1 = L
               L2 = L
            } else {
               if ( B( L, L-1 ).NE.ZERO ) {
                  L1 = L - 1
                  L2 = L
                  LNEXT = L - 2
               } else {
                  L1 = L
                  L2 = L
                  LNEXT = L - 1
               }
            }

            // Start row loop (index = K)
            // K1 (K2): row index of the first (last) row of X(K,L)

            KNEXT = 1
            DO 170 K = 1, M
               IF( K.LT.KNEXT ) GO TO 170
               if ( K.EQ.M ) {
                  K1 = K
                  K2 = K
               } else {
                  if ( A( K+1, K ).NE.ZERO ) {
                     K1 = K
                     K2 = K + 1
                     KNEXT = K + 2
                  } else {
                     K1 = K
                     K2 = K
                     KNEXT = K + 1
                  }
               }

               if ( L1.EQ.L2 .AND. K1.EQ.K2 ) {
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L1, C( K1, MIN( L1+1, N ) ), LDC, B( L1, MIN( L1+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE

                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  if ( DA11.LE.SMIN ) {
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  }
                  DB = ABS( VEC( 1, 1 ) )
                  if ( DA11.LT.ONE .AND. DB.GT.ONE ) {
                     IF( DB.GT.BIGNUM*DA11 ) SCALOC = ONE / DB
                  }
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11

                  if ( SCALOC.NE.ONE ) {
                     DO 130 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
  130                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )

               } else if ( L1.EQ.L2 .AND. K1.NE.K2 ) {

                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )

                  dlaln2(.TRUE., 2, 1, SMIN, ONE, A( K1, K1 ), LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  IF( IERR.NE.0 ) INFO = 1

                  if ( SCALOC.NE.ONE ) {
                     DO 140 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
  140                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )

               } else if ( L1.NE.L2 .AND. K1.EQ.K2 ) {

                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )

                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )

                  dlaln2(.FALSE., 2, 1, SMIN, ONE, B( L1, L1 ), LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  IF( IERR.NE.0 ) INFO = 1

                  if ( SCALOC.NE.ONE ) {
                     DO 150 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
  150                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )

               } else if ( L1.NE.L2 .AND. K1.NE.K2 ) {

                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )

                  dlasy2(.TRUE., .TRUE., ISGN, 2, 2, A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, 2, XNORM, IERR );
                  IF( IERR.NE.0 ) INFO = 1

                  if ( SCALOC.NE.ONE ) {
                     DO 160 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
  160                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               }

  170       CONTINUE
  180    CONTINUE

      } else if ( NOTRNA .AND. .NOT.NOTRNB ) {

         // Solve    A*X + ISGN*X*B**T = scale*C.

         // The (K,L)th block of X is determined starting from
         // bottom-right corner column by column by

             // A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)

         // Where
                       // M                          N
             // R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
                     // I=K+1                      J=L+1

         // Start column loop (index = L)
         // L1 (L2): column index of the first (last) row of X(K,L)

         LNEXT = N
         DO 240 L = N, 1, -1
            IF( L.GT.LNEXT ) GO TO 240
            if ( L.EQ.1 ) {
               L1 = L
               L2 = L
            } else {
               if ( B( L, L-1 ).NE.ZERO ) {
                  L1 = L - 1
                  L2 = L
                  LNEXT = L - 2
               } else {
                  L1 = L
                  L2 = L
                  LNEXT = L - 1
               }
            }

            // Start row loop (index = K)
            // K1 (K2): row index of the first (last) row of X(K,L)

            KNEXT = M
            DO 230 K = M, 1, -1
               IF( K.GT.KNEXT ) GO TO 230
               if ( K.EQ.1 ) {
                  K1 = K
                  K2 = K
               } else {
                  if ( A( K, K-1 ).NE.ZERO ) {
                     K1 = K - 1
                     K2 = K
                     KNEXT = K - 2
                  } else {
                     K1 = K
                     K2 = K
                     KNEXT = K - 1
                  }
               }

               if ( L1.EQ.L2 .AND. K1.EQ.K2 ) {
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, C( MIN( K1+1, M ), L1 ), 1 )                   SUMR = DDOT( N-L1, C( K1, MIN( L1+1, N ) ), LDC, B( L1, MIN( L1+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE

                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  if ( DA11.LE.SMIN ) {
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  }
                  DB = ABS( VEC( 1, 1 ) )
                  if ( DA11.LT.ONE .AND. DB.GT.ONE ) {
                     IF( DB.GT.BIGNUM*DA11 ) SCALOC = ONE / DB
                  }
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11

                  if ( SCALOC.NE.ONE ) {
                     DO 190 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
  190                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )

               } else if ( L1.EQ.L2 .AND. K1.NE.K2 ) {

                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, C( MIN( K2+1, M ), L1 ), 1 )                   SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, C( MIN( K2+1, M ), L1 ), 1 )                   SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )

                  dlaln2(.FALSE., 2, 1, SMIN, ONE, A( K1, K1 ), LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  IF( IERR.NE.0 ) INFO = 1

                  if ( SCALOC.NE.ONE ) {
                     DO 200 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
  200                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )

               } else if ( L1.NE.L2 .AND. K1.EQ.K2 ) {

                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, C( MIN( K1+1, M ), L1 ), 1 )                   SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )

                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, C( MIN( K1+1, M ), L2 ), 1 )                   SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )

                  dlaln2(.FALSE., 2, 1, SMIN, ONE, B( L1, L1 ), LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  IF( IERR.NE.0 ) INFO = 1

                  if ( SCALOC.NE.ONE ) {
                     DO 210 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
  210                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )

               } else if ( L1.NE.L2 .AND. K1.NE.K2 ) {

                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, C( MIN( K2+1, M ), L1 ), 1 )                   SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, C( MIN( K2+1, M ), L2 ), 1 )                   SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, C( MIN( K2+1, M ), L1 ), 1 )                   SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )

                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, C( MIN( K2+1, M ), L2 ), 1 )                   SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )

                  dlasy2(.FALSE., .TRUE., ISGN, 2, 2, A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, 2, XNORM, IERR );
                  IF( IERR.NE.0 ) INFO = 1

                  if ( SCALOC.NE.ONE ) {
                     DO 220 J = 1, N
                        dscal(M, SCALOC, C( 1, J ), 1 );
  220                CONTINUE
                     SCALE = SCALE*SCALOC
                  }
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               }

  230       CONTINUE
  240    CONTINUE

      }

      RETURN

      // End of DTRSYL

      }
