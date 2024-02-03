      SUBROUTINE CTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, LDC, SCALE, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANA, TRANB;
      int                INFO, ISGN, LDA, LDB, LDC, M, N;
      REAL               SCALE
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRNA, NOTRNB;
      int                J, K, L;
      REAL               BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN, SMLNUM
      COMPLEX            A11, SUML, SUMR, VEC, X11
      // ..
      // .. Local Arrays ..
      REAL               DUM( 1 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, SLAMCH
      COMPLEX            CDOTC, CDOTU, CLADIV
      // EXTERNAL LSAME, CLANGE, SLAMCH, CDOTC, CDOTU, CLADIV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, CONJG, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Decode and Test input parameters

      NOTRNA = LSAME( TRANA, 'N' )
      NOTRNB = LSAME( TRANB, 'N' )

      INFO = 0
      if ( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'C' ) ) {
         INFO = -1
      } else if ( .NOT.NOTRNB .AND. .NOT.LSAME( TRANB, 'C' ) ) {
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
         xerbla('CTRSYL', -INFO );
         RETURN
      }

      // Quick return if possible

      SCALE = ONE
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      // Set constants to control overflow

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SMLNUM*REAL( M*N ) / EPS
      BIGNUM = ONE / SMLNUM
      SMIN = MAX( SMLNUM, EPS*CLANGE( 'M', M, M, A, LDA, DUM ), EPS*CLANGE( 'M', N, N, B, LDB, DUM ) )
      SGN = ISGN

      if ( NOTRNA .AND. NOTRNB ) {

         // Solve    A*X + ISGN*X*B = scale*C.

         // The (K,L)th block of X is determined starting from
         // bottom-left corner column by column by

             // A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

         // Where
                     // M                        L-1
           // R(K,L) = SUM [A(K,I)*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)].
                   // I=K+1                      J=1

         DO 30 L = 1, N
            DO 20 K = M, 1, -1

               SUML = CDOTU( M-K, A( K, MIN( K+1, M ) ), LDA, C( MIN( K+1, M ), L ), 1 )
               SUMR = CDOTU( L-1, C( K, 1 ), LDC, B( 1, L ), 1 )
               VEC = C( K, L ) - ( SUML+SGN*SUMR )

               SCALOC = ONE
               A11 = A( K, K ) + SGN*B( L, L )
               DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
               if ( DA11.LE.SMIN ) {
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               }
               DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
               if ( DA11.LT.ONE .AND. DB.GT.ONE ) {
                  IF( DB.GT.BIGNUM*DA11 ) SCALOC = ONE / DB
               }
               X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )

               if ( SCALOC.NE.ONE ) {
                  DO 10 J = 1, N
                     csscal(M, SCALOC, C( 1, J ), 1 );
   10             CONTINUE
                  SCALE = SCALE*SCALOC
               }
               C( K, L ) = X11

   20       CONTINUE
   30    CONTINUE

      } else if ( .NOT.NOTRNA .AND. NOTRNB ) {

         // Solve    A**H *X + ISGN*X*B = scale*C.

         // The (K,L)th block of X is determined starting from
         // upper-left corner column by column by

             // A**H(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

         // Where
                    // K-1                           L-1
           // R(K,L) = SUM [A**H(I,K)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]
                    // I=1                           J=1

         DO 60 L = 1, N
            DO 50 K = 1, M

               SUML = CDOTC( K-1, A( 1, K ), 1, C( 1, L ), 1 )
               SUMR = CDOTU( L-1, C( K, 1 ), LDC, B( 1, L ), 1 )
               VEC = C( K, L ) - ( SUML+SGN*SUMR )

               SCALOC = ONE
               A11 = CONJG( A( K, K ) ) + SGN*B( L, L )
               DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
               if ( DA11.LE.SMIN ) {
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               }
               DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
               if ( DA11.LT.ONE .AND. DB.GT.ONE ) {
                  IF( DB.GT.BIGNUM*DA11 ) SCALOC = ONE / DB
               }

               X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )

               if ( SCALOC.NE.ONE ) {
                  DO 40 J = 1, N
                     csscal(M, SCALOC, C( 1, J ), 1 );
   40             CONTINUE
                  SCALE = SCALE*SCALOC
               }
               C( K, L ) = X11

   50       CONTINUE
   60    CONTINUE

      } else if ( .NOT.NOTRNA .AND. .NOT.NOTRNB ) {

         // Solve    A**H*X + ISGN*X*B**H = C.

         // The (K,L)th block of X is determined starting from
         // upper-right corner column by column by

             // A**H(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)

         // Where
                     // K-1
            // R(K,L) = SUM [A**H(I,K)*X(I,L)] +
                     // I=1
                            // N
                      // ISGN*SUM [X(K,J)*B**H(L,J)].
                           // J=L+1

         DO 90 L = N, 1, -1
            DO 80 K = 1, M

               SUML = CDOTC( K-1, A( 1, K ), 1, C( 1, L ), 1 )
               SUMR = CDOTC( N-L, C( K, MIN( L+1, N ) ), LDC, B( L, MIN( L+1, N ) ), LDB )
               VEC = C( K, L ) - ( SUML+SGN*CONJG( SUMR ) )

               SCALOC = ONE
               A11 = CONJG( A( K, K )+SGN*B( L, L ) )
               DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
               if ( DA11.LE.SMIN ) {
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               }
               DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
               if ( DA11.LT.ONE .AND. DB.GT.ONE ) {
                  IF( DB.GT.BIGNUM*DA11 ) SCALOC = ONE / DB
               }

               X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )

               if ( SCALOC.NE.ONE ) {
                  DO 70 J = 1, N
                     csscal(M, SCALOC, C( 1, J ), 1 );
   70             CONTINUE
                  SCALE = SCALE*SCALOC
               }
               C( K, L ) = X11

   80       CONTINUE
   90    CONTINUE

      } else if ( NOTRNA .AND. .NOT.NOTRNB ) {

         // Solve    A*X + ISGN*X*B**H = C.

         // The (K,L)th block of X is determined starting from
         // bottom-left corner column by column by

            // A(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)

         // Where
                     // M                          N
           // R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B**H(L,J)]
                   // I=K+1                      J=L+1

         DO 120 L = N, 1, -1
            DO 110 K = M, 1, -1

               SUML = CDOTU( M-K, A( K, MIN( K+1, M ) ), LDA, C( MIN( K+1, M ), L ), 1 )                SUMR = CDOTC( N-L, C( K, MIN( L+1, N ) ), LDC, B( L, MIN( L+1, N ) ), LDB )
               VEC = C( K, L ) - ( SUML+SGN*CONJG( SUMR ) )

               SCALOC = ONE
               A11 = A( K, K ) + SGN*CONJG( B( L, L ) )
               DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
               if ( DA11.LE.SMIN ) {
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               }
               DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
               if ( DA11.LT.ONE .AND. DB.GT.ONE ) {
                  IF( DB.GT.BIGNUM*DA11 ) SCALOC = ONE / DB
               }

               X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )

               if ( SCALOC.NE.ONE ) {
                  DO 100 J = 1, N
                     csscal(M, SCALOC, C( 1, J ), 1 );
  100             CONTINUE
                  SCALE = SCALE*SCALOC
               }
               C( K, L ) = X11

  110       CONTINUE
  120    CONTINUE

      }

      RETURN

      // End of CTRSYL

      }
