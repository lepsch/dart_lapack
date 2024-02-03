      void dbdsdc(UPLO, COMPQ, N, D, E, U, LDU, VT, LDVT, Q, IQ, WORK, IWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPQ, UPLO;
      int                INFO, LDU, LDVT, N;
      // ..
      // .. Array Arguments ..
      int                IQ( * ), IWORK( * );
      double             D( * ), E( * ), Q( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

// =====================================================================
// Changed dimension statement in comment describing E from (N) to
// (N-1).  Sven, 17 Feb 05.
// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      int                DIFL, DIFR, GIVCOL, GIVNUM, GIVPTR, I, IC, ICOMPQ, IERR, II, IS, IU, IUPLO, IVT, J, K, KK, MLVL, NM1, NSIZE, PERM, POLES, QSTART, SMLSIZ, SMLSZP, SQRE, START, WSTART, Z;
      double             CS, EPS, ORGNRM, P, R, SN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, DLANST;
      // EXTERNAL LSAME, ILAENV, DLAMCH, DLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLARTG, DLASCL, DLASD0, DLASDA, DLASDQ, DLASET, DLASR, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, INT, LOG, SIGN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;

      IUPLO = 0;
      if( LSAME( UPLO, 'U' ) ) IUPLO = 1;
      IF( LSAME( UPLO, 'L' ) ) IUPLO = 2;
      if ( LSAME( COMPQ, 'N' ) ) {
         ICOMPQ = 0;
      } else if ( LSAME( COMPQ, 'P' ) ) {
         ICOMPQ = 1;
      } else if ( LSAME( COMPQ, 'I' ) ) {
         ICOMPQ = 2;
      } else {
         ICOMPQ = -1;
      }
      if ( IUPLO == 0 ) {
         INFO = -1;
      } else if ( ICOMPQ < 0 ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( ( LDU < 1 ) || ( ( ICOMPQ == 2 ) && ( LDU < N ) ) ) {
         INFO = -7;
      } else if ( ( LDVT < 1 ) || ( ( ICOMPQ == 2 ) && ( LDVT < N ) ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('DBDSDC', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;
      SMLSIZ = ILAENV( 9, 'DBDSDC', ' ', 0, 0, 0, 0 );
      if ( N == 1 ) {
         if ( ICOMPQ == 1 ) {
            Q( 1 ) = SIGN( ONE, D( 1 ) );
            Q( 1+SMLSIZ*N ) = ONE;
         } else if ( ICOMPQ == 2 ) {
            U( 1, 1 ) = SIGN( ONE, D( 1 ) );
            VT( 1, 1 ) = ONE;
         }
         D( 1 ) = ABS( D( 1 ) );
         return;
      }
      NM1 = N - 1;

      // If matrix lower bidiagonal, rotate to be upper bidiagonal
      // by applying Givens rotations on the left

      WSTART = 1;
      QSTART = 3;
      if ( ICOMPQ == 1 ) {
         dcopy(N, D, 1, Q( 1 ), 1 );
         dcopy(N-1, E, 1, Q( N+1 ), 1 );
      }
      if ( IUPLO == 2 ) {
         QSTART = 5;
         if (ICOMPQ == 2) WSTART = 2*N - 1;
         for (I = 1; I <= N - 1; I++) { // 10
            dlartg(D( I ), E( I ), CS, SN, R );
            D( I ) = R;
            E( I ) = SN*D( I+1 );
            D( I+1 ) = CS*D( I+1 );
            if ( ICOMPQ == 1 ) {
               Q( I+2*N ) = CS;
               Q( I+3*N ) = SN;
            } else if ( ICOMPQ == 2 ) {
               WORK( I ) = CS;
               WORK( NM1+I ) = -SN;
            }
         } // 10
      }

      // If ICOMPQ = 0, use DLASDQ to compute the singular values.

      if ( ICOMPQ == 0 ) {
         // Ignore WSTART, instead using WORK( 1 ), since the two vectors
         // for CS and -SN above are added only if ICOMPQ == 2,
         // and adding them exceeds documented WORK size of 4*n.
         dlasdq('U', 0, N, 0, 0, 0, D, E, VT, LDVT, U, LDU, U, LDU, WORK( 1 ), INFO );
         GO TO 40;
      }

      // If N is smaller than the minimum divide size SMLSIZ, then solve
      // the problem with another solver.

      if ( N <= SMLSIZ ) {
         if ( ICOMPQ == 2 ) {
            dlaset('A', N, N, ZERO, ONE, U, LDU );
            dlaset('A', N, N, ZERO, ONE, VT, LDVT );
            dlasdq('U', 0, N, N, N, 0, D, E, VT, LDVT, U, LDU, U, LDU, WORK( WSTART ), INFO );
         } else if ( ICOMPQ == 1 ) {
            IU = 1;
            IVT = IU + N;
            dlaset('A', N, N, ZERO, ONE, Q( IU+( QSTART-1 )*N ), N );
            dlaset('A', N, N, ZERO, ONE, Q( IVT+( QSTART-1 )*N ), N );
            dlasdq('U', 0, N, N, N, 0, D, E, Q( IVT+( QSTART-1 )*N ), N, Q( IU+( QSTART-1 )*N ), N, Q( IU+( QSTART-1 )*N ), N, WORK( WSTART ), INFO );
         }
         GO TO 40;
      }

      if ( ICOMPQ == 2 ) {
         dlaset('A', N, N, ZERO, ONE, U, LDU );
         dlaset('A', N, N, ZERO, ONE, VT, LDVT );
      }

      // Scale.

      ORGNRM = DLANST( 'M', N, D, E );
      if (ORGNRM == ZERO) return;
      dlascl('G', 0, 0, ORGNRM, ONE, N, 1, D, N, IERR );
      dlascl('G', 0, 0, ORGNRM, ONE, NM1, 1, E, NM1, IERR );

      EPS = (0.9)*DLAMCH( 'Epsilon' );

      MLVL = INT( LOG( DBLE( N ) / DBLE( SMLSIZ+1 ) ) / LOG( TWO ) ) + 1;
      SMLSZP = SMLSIZ + 1;

      if ( ICOMPQ == 1 ) {
         IU = 1;
         IVT = 1 + SMLSIZ;
         DIFL = IVT + SMLSZP;
         DIFR = DIFL + MLVL;
         Z = DIFR + MLVL*2;
         IC = Z + MLVL;
         IS = IC + 1;
         POLES = IS + 1;
         GIVNUM = POLES + 2*MLVL;

         K = 1;
         GIVPTR = 2;
         PERM = 3;
         GIVCOL = PERM + MLVL;
      }

      for (I = 1; I <= N; I++) { // 20
         if ( ABS( D( I ) ) < EPS ) {
            D( I ) = SIGN( EPS, D( I ) );
         }
      } // 20

      START = 1;
      SQRE = 0;

      for (I = 1; I <= NM1; I++) { // 30
         if ( ( ABS( E( I ) ) < EPS ) || ( I == NM1 ) ) {

            // Subproblem found. First determine its size and then
            // apply divide and conquer on it.

            if ( I < NM1 ) {

               // A subproblem with E(I) small for I < NM1.

               NSIZE = I - START + 1;
            } else if ( ABS( E( I ) ) >= EPS ) {

               // A subproblem with E(NM1) not too small but I = NM1.

               NSIZE = N - START + 1;
            } else {

               // A subproblem with E(NM1) small. This implies an
               // 1-by-1 subproblem at D(N). Solve this 1-by-1 problem
               // first.

               NSIZE = I - START + 1;
               if ( ICOMPQ == 2 ) {
                  U( N, N ) = SIGN( ONE, D( N ) );
                  VT( N, N ) = ONE;
               } else if ( ICOMPQ == 1 ) {
                  Q( N+( QSTART-1 )*N ) = SIGN( ONE, D( N ) );
                  Q( N+( SMLSIZ+QSTART-1 )*N ) = ONE;
               }
               D( N ) = ABS( D( N ) );
            }
            if ( ICOMPQ == 2 ) {
               dlasd0(NSIZE, SQRE, D( START ), E( START ), U( START, START ), LDU, VT( START, START ), LDVT, SMLSIZ, IWORK, WORK( WSTART ), INFO );
            } else {
               dlasda(ICOMPQ, SMLSIZ, NSIZE, SQRE, D( START ), E( START ), Q( START+( IU+QSTART-2 )*N ), N, Q( START+( IVT+QSTART-2 )*N ), IQ( START+K*N ), Q( START+( DIFL+QSTART-2 )* N ), Q( START+( DIFR+QSTART-2 )*N ), Q( START+( Z+QSTART-2 )*N ), Q( START+( POLES+QSTART-2 )*N ), IQ( START+GIVPTR*N ), IQ( START+GIVCOL*N ), N, IQ( START+PERM*N ), Q( START+( GIVNUM+QSTART-2 )*N ), Q( START+( IC+QSTART-2 )*N ), Q( START+( IS+QSTART-2 )*N ), WORK( WSTART ), IWORK, INFO );
            }
            if ( INFO != 0 ) {
               return;
            }
            START = I + 1;
         }
      } // 30

      // Unscale

      dlascl('G', 0, 0, ONE, ORGNRM, N, 1, D, N, IERR );
      } // 40

      // Use Selection Sort to minimize swaps of singular vectors

      for (II = 2; II <= N; II++) { // 60
         I = II - 1;
         KK = I;
         P = D( I );
         for (J = II; J <= N; J++) { // 50
            if ( D( J ) > P ) {
               KK = J;
               P = D( J );
            }
         } // 50
         if ( KK != I ) {
            D( KK ) = D( I );
            D( I ) = P;
            if ( ICOMPQ == 1 ) {
               IQ( I ) = KK;
            } else if ( ICOMPQ == 2 ) {
               dswap(N, U( 1, I ), 1, U( 1, KK ), 1 );
               dswap(N, VT( I, 1 ), LDVT, VT( KK, 1 ), LDVT );
            }
         } else if ( ICOMPQ == 1 ) {
            IQ( I ) = I;
         }
      } // 60

      // If ICOMPQ = 1, use IQ(N,1) as the indicator for UPLO

      if ( ICOMPQ == 1 ) {
         if ( IUPLO == 1 ) {
            IQ( N ) = 1;
         } else {
            IQ( N ) = 0;
         }
      }

      // If B is lower bidiagonal, update U by those Givens rotations
      // which rotated B to be upper bidiagonal

      if( ( IUPLO == 2 ) && ( ICOMPQ == 2 ) ) dlasr( 'L', 'V', 'B', N, N, WORK( 1 ), WORK( N ), U, LDU );

      return;
      }
