      void zhgeqz(JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT, ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, RWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPQ, COMPZ, JOB;
      int                IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      Complex         ALPHA( * ), BETA( * ), H( LDH, * ), Q( LDQ, * ), T( LDT, * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             HALF;
      const              HALF = 0.5 ;
      // ..
      // .. Local Scalars ..
      bool               ILAZR2, ILAZRO, ILQ, ILSCHR, ILZ, LQUERY;
      int                ICOMPQ, ICOMPZ, IFIRST, IFRSTM, IITER, ILAST, ILASTM, IN, ISCHUR, ISTART, J, JC, JCH, JITER, JR, MAXIT;
      double             ABSB, ANORM, ASCALE, ATOL, BNORM, BSCALE, BTOL, C, SAFMIN, TEMP, TEMP2, TEMPR, ULP;
      Complex         ABI22, AD11, AD12, AD21, AD22, CTEMP, CTEMP2, CTEMP3, ESHIFT, S, SHIFT, SIGNBC, U12, X, ABI12, Y;
      // ..
      // .. External Functions ..
      //- Complex         ZLADIV;
      //- bool               lsame;
      //- double             DLAMCH, ZLANHS;
      // EXTERNAL ZLADIV, lsame, DLAMCH, ZLANHS
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLARTG, ZLASET, ZROT, ZSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN, SQRT
      // ..
      // .. Statement Functions ..
      double             ABS1;
      // ..
      // .. Statement Function definitions ..
      ABS1[X] = ( X.toDouble() ).abs() + ( DIMAG( X ) ).abs();
      // ..
      // .. Executable Statements ..

      // Decode JOB, COMPQ, COMPZ

      if ( lsame( JOB, 'E' ) ) {
         ILSCHR = false;
         ISCHUR = 1;
      } else if ( lsame( JOB, 'S' ) ) {
         ILSCHR = true;
         ISCHUR = 2;
      } else {
         ILSCHR = true;
         ISCHUR = 0;
      }

      if ( lsame( COMPQ, 'N' ) ) {
         ILQ = false;
         ICOMPQ = 1;
      } else if ( lsame( COMPQ, 'V' ) ) {
         ILQ = true;
         ICOMPQ = 2;
      } else if ( lsame( COMPQ, 'I' ) ) {
         ILQ = true;
         ICOMPQ = 3;
      } else {
         ILQ = true;
         ICOMPQ = 0;
      }

      if ( lsame( COMPZ, 'N' ) ) {
         ILZ = false;
         ICOMPZ = 1;
      } else if ( lsame( COMPZ, 'V' ) ) {
         ILZ = true;
         ICOMPZ = 2;
      } else if ( lsame( COMPZ, 'I' ) ) {
         ILZ = true;
         ICOMPZ = 3;
      } else {
         ILZ = true;
         ICOMPZ = 0;
      }

      // Check Argument Values

      INFO = 0;
      WORK[1] = max( 1, N );
      LQUERY = ( LWORK == -1 );
      if ( ISCHUR == 0 ) {
         INFO = -1;
      } else if ( ICOMPQ == 0 ) {
         INFO = -2;
      } else if ( ICOMPZ == 0 ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( ILO < 1 ) {
         INFO = -5;
      } else if ( IHI > N || IHI < ILO-1 ) {
         INFO = -6;
      } else if ( LDH < N ) {
         INFO = -8;
      } else if ( LDT < N ) {
         INFO = -10;
      } else if ( LDQ < 1 || ( ILQ && LDQ < N ) ) {
         INFO = -14;
      } else if ( LDZ < 1 || ( ILZ && LDZ < N ) ) {
         INFO = -16;
      } else if ( LWORK < max( 1, N ) && !LQUERY ) {
         INFO = -18;
      }
      if ( INFO != 0 ) {
         xerbla('ZHGEQZ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      // WORK( 1 ) = CMPLX( 1 )
      if ( N <= 0 ) {
         WORK[1] = DCMPLX( 1 );
         return;
      }

      // Initialize Q and Z

      if (ICOMPQ == 3) zlaset( 'Full', N, N, CZERO, CONE, Q, LDQ );
      IF( ICOMPZ == 3 ) zlaset( 'Full', N, N, CZERO, CONE, Z, LDZ );

      // Machine Constants

      IN = IHI + 1 - ILO;
      SAFMIN = DLAMCH( 'S' );
      ULP = DLAMCH( 'E' )*DLAMCH( 'B' );
      ANORM = ZLANHS( 'F', IN, H( ILO, ILO ), LDH, RWORK );
      BNORM = ZLANHS( 'F', IN, T( ILO, ILO ), LDT, RWORK );
      ATOL = max( SAFMIN, ULP*ANORM );
      BTOL = max( SAFMIN, ULP*BNORM );
      ASCALE = ONE / max( SAFMIN, ANORM );
      BSCALE = ONE / max( SAFMIN, BNORM );


      // Set Eigenvalues IHI+1:N

      for (J = IHI + 1; J <= N; J++) { // 10
         ABSB = ( T( J, J ) ).abs();
         if ( ABSB > SAFMIN ) {
            SIGNBC = DCONJG( T( J, J ) / ABSB );
            T[J, J] = ABSB;
            if ( ILSCHR ) {
               zscal(J-1, SIGNBC, T( 1, J ), 1 );
               zscal(J, SIGNBC, H( 1, J ), 1 );
            } else {
               zscal(1, SIGNBC, H( J, J ), 1 );
            }
            if (ILZ) zscal( N, SIGNBC, Z( 1, J ), 1 );
         } else {
            T[J, J] = CZERO;
         }
         ALPHA[J] = H( J, J );
         BETA[J] = T( J, J );
      } // 10

      // If IHI < ILO, skip QZ steps

      if (IHI < ILO) GO TO 190;

      // MAIN QZ ITERATION LOOP

      // Initialize dynamic indices

      // Eigenvalues ILAST+1:N have been found.
         // Column operations modify rows IFRSTM:whatever
         // Row operations modify columns whatever:ILASTM

      // If only eigenvalues are being computed, then
         // IFRSTM is the row of the last splitting row above row ILAST;
         // this is always at least ILO.
      // IITER counts iterations since the last eigenvalue was found,
         // to tell when to use an extraordinary shift.
      // MAXIT is the maximum number of QZ sweeps allowed.

      ILAST = IHI;
      if ( ILSCHR ) {
         IFRSTM = 1;
         ILASTM = N;
      } else {
         IFRSTM = ILO;
         ILASTM = IHI;
      }
      IITER = 0;
      ESHIFT = CZERO;
      MAXIT = 30*( IHI-ILO+1 );

      for (JITER = 1; JITER <= MAXIT; JITER++) { // 170

         // Check for too many iterations.

         if (JITER > MAXIT) GO TO 180;

         // Split the matrix if possible.

         // Two tests:
            // 1: H(j,j-1)=0  or  j=ILO
            // 2: T(j,j)=0

         // Special case: j=ILAST

         if ( ILAST == ILO ) {
            GO TO 60;
         } else {
            if ( ABS1( H( ILAST, ILAST-1 ) ) <= max( SAFMIN, ULP*(  ABS1( H( ILAST, ILAST ) ) + ABS1( H( ILAST-1, ILAST-1 ) ) ) ) ) {
               H[ILAST, ILAST-1] = CZERO;
               GO TO 60;
            }
         }

         if ( ( T( ILAST, ILAST ) ).abs() <= BTOL ) {
            T[ILAST, ILAST] = CZERO;
            GO TO 50;
         }

         // General case: j<ILAST

         for (J = ILAST - 1; J >= ILO; J--) { // 40

            // Test 1: for H(j,j-1)=0 or j=ILO

            if ( J == ILO ) {
               ILAZRO = true;
            } else {
               if ( ABS1( H( J, J-1 ) ) <= max( SAFMIN, ULP*(  ABS1( H( J, J ) ) + ABS1( H( J-1, J-1 ) ) ) ) ) {
                  H[J, J-1] = CZERO;
                  ILAZRO = true;
               } else {
                  ILAZRO = false;
               }
            }

            // Test 2: for T(j,j)=0

            if ( ( T( J, J ) ).abs() < BTOL ) {
               T[J, J] = CZERO;

               // Test 1a: Check for 2 consecutive small subdiagonals in A

               ILAZR2 = false;
               if ( !ILAZRO ) {
                  if( ABS1( H( J, J-1 ) )*( ASCALE*ABS1( H( J+1, J ) ) ) <= ABS1( H( J, J ) )*( ASCALE*ATOL ) ) ILAZR2 = true;
               }

               // If both tests pass (1 & 2), i.e., the leading diagonal
               // element of B in the block is zero, split a 1x1 block off
               // at the top. (I.e., at the J-th row/column) The leading
               // diagonal element of the remainder can also be zero, so
               // this may have to be done repeatedly.

               if ( ILAZRO || ILAZR2 ) {
                  for (JCH = J; JCH <= ILAST - 1; JCH++) { // 20
                     CTEMP = H( JCH, JCH );
                     zlartg(CTEMP, H( JCH+1, JCH ), C, S, H( JCH, JCH ) );
                     H[JCH+1, JCH] = CZERO;
                     zrot(ILASTM-JCH, H( JCH, JCH+1 ), LDH, H( JCH+1, JCH+1 ), LDH, C, S );
                     zrot(ILASTM-JCH, T( JCH, JCH+1 ), LDT, T( JCH+1, JCH+1 ), LDT, C, S )                      IF( ILQ ) CALL ZROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1, C, DCONJG( S ) );
                     if (ILAZR2) H( JCH, JCH-1 ) = H( JCH, JCH-1 )*C;
                     ILAZR2 = false;
                     if ( ABS1( T( JCH+1, JCH+1 ) ) >= BTOL ) {
                        if ( JCH+1 >= ILAST ) {
                           GO TO 60;
                        } else {
                           IFIRST = JCH + 1;
                           GO TO 70;
                        }
                     }
                     T[JCH+1, JCH+1] = CZERO;
                  } // 20
                  GO TO 50;
               } else {

                  // Only test 2 passed -- chase the zero to T(ILAST,ILAST)
                  // Then process as in the case T(ILAST,ILAST)=0

                  for (JCH = J; JCH <= ILAST - 1; JCH++) { // 30
                     CTEMP = T( JCH, JCH+1 );
                     zlartg(CTEMP, T( JCH+1, JCH+1 ), C, S, T( JCH, JCH+1 ) );
                     T[JCH+1, JCH+1] = CZERO;
                     if (JCH < ILASTM-1) zrot( ILASTM-JCH-1, T( JCH, JCH+2 ), LDT, T( JCH+1, JCH+2 ), LDT, C, S );
                     zrot(ILASTM-JCH+2, H( JCH, JCH-1 ), LDH, H( JCH+1, JCH-1 ), LDH, C, S )                      IF( ILQ ) CALL ZROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1, C, DCONJG( S ) );
                     CTEMP = H( JCH+1, JCH );
                     zlartg(CTEMP, H( JCH+1, JCH-1 ), C, S, H( JCH+1, JCH ) );
                     H[JCH+1, JCH-1] = CZERO;
                     zrot(JCH+1-IFRSTM, H( IFRSTM, JCH ), 1, H( IFRSTM, JCH-1 ), 1, C, S );
                     zrot(JCH-IFRSTM, T( IFRSTM, JCH ), 1, T( IFRSTM, JCH-1 ), 1, C, S )                      IF( ILZ ) CALL ZROT( N, Z( 1, JCH ), 1, Z( 1, JCH-1 ), 1, C, S );
                  } // 30
                  GO TO 50;
               }
            } else if ( ILAZRO ) {

               // Only test 1 passed -- work on J:ILAST

               IFIRST = J;
               GO TO 70;
            }

            // Neither test passed -- try next J

         } // 40

         // (Drop-through is "impossible")

         INFO = 2*N + 1;
         GO TO 210;

         // T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
         // 1x1 block.

         } // 50
         CTEMP = H( ILAST, ILAST );
         zlartg(CTEMP, H( ILAST, ILAST-1 ), C, S, H( ILAST, ILAST ) );
         H[ILAST, ILAST-1] = CZERO;
         zrot(ILAST-IFRSTM, H( IFRSTM, ILAST ), 1, H( IFRSTM, ILAST-1 ), 1, C, S );
         zrot(ILAST-IFRSTM, T( IFRSTM, ILAST ), 1, T( IFRSTM, ILAST-1 ), 1, C, S )          IF( ILZ ) CALL ZROT( N, Z( 1, ILAST ), 1, Z( 1, ILAST-1 ), 1, C, S );

         // H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHA and BETA

         } // 60
         ABSB = ( T( ILAST, ILAST ) ).abs();
         if ( ABSB > SAFMIN ) {
            SIGNBC = DCONJG( T( ILAST, ILAST ) / ABSB );
            T[ILAST, ILAST] = ABSB;
            if ( ILSCHR ) {
               zscal(ILAST-IFRSTM, SIGNBC, T( IFRSTM, ILAST ), 1 );
               zscal(ILAST+1-IFRSTM, SIGNBC, H( IFRSTM, ILAST ), 1 );
            } else {
               zscal(1, SIGNBC, H( ILAST, ILAST ), 1 );
            }
            if (ILZ) zscal( N, SIGNBC, Z( 1, ILAST ), 1 );
         } else {
            T[ILAST, ILAST] = CZERO;
         }
         ALPHA[ILAST] = H( ILAST, ILAST );
         BETA[ILAST] = T( ILAST, ILAST );

         // Go to next block -- exit if finished.

         ILAST = ILAST - 1;
         if (ILAST < ILO) GO TO 190;

         // Reset counters

         IITER = 0;
         ESHIFT = CZERO;
         if ( !ILSCHR ) {
            ILASTM = ILAST;
            if (IFRSTM > ILAST) IFRSTM = ILO;
         }
         GO TO 160;

         // QZ step

         // This iteration only involves rows/columns IFIRST:ILAST.  We
         // assume IFIRST < ILAST, and that the diagonal of B is non-zero.

         } // 70
         IITER = IITER + 1;
         if ( !ILSCHR ) {
            IFRSTM = IFIRST;
         }

         // Compute the Shift.

         // At this point, IFIRST < ILAST, and the diagonal elements of
         // T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
         // magnitude)

         if ( ( IITER / 10 )*10 != IITER ) {

            // The Wilkinson shift (AEP p.512), i.e., the eigenvalue of
            // the bottom-right 2x2 block of A inv(B) which is nearest to
            // the bottom-right element.

            // We factor B as U*D, where U has unit diagonals, and
            // compute (A*inv(D))*inv(U).

            U12 = ( BSCALE*T( ILAST-1, ILAST ) ) / ( BSCALE*T( ILAST, ILAST ) )             AD11 = ( ASCALE*H( ILAST-1, ILAST-1 ) ) / ( BSCALE*T( ILAST-1, ILAST-1 ) )             AD21 = ( ASCALE*H( ILAST, ILAST-1 ) ) / ( BSCALE*T( ILAST-1, ILAST-1 ) )             AD12 = ( ASCALE*H( ILAST-1, ILAST ) ) / ( BSCALE*T( ILAST, ILAST ) )             AD22 = ( ASCALE*H( ILAST, ILAST ) ) / ( BSCALE*T( ILAST, ILAST ) );
            ABI22 = AD22 - U12*AD21;
            ABI12 = AD12 - U12*AD11;

            SHIFT = ABI22;
            CTEMP = sqrt( ABI12 )*sqrt( AD21 );
            TEMP = ABS1( CTEMP );
            if ( CTEMP != ZERO ) {
               X = HALF*( AD11-SHIFT );
               TEMP2 = ABS1( X );
               TEMP = max( TEMP, ABS1( X ) );
               Y = TEMP*sqrt( ( X / TEMP )**2+( CTEMP / TEMP )**2 );
               if ( TEMP2 > ZERO ) {
                  if( (X / TEMP2).toDouble()*Y.toDouble()+ DIMAG( X / TEMP2 )*DIMAG( Y ) < ZERO )Y = -Y;
               }
               SHIFT = SHIFT - CTEMP*ZLADIV( CTEMP, ( X+Y ) );
            }
         } else {

            // Exceptional shift.  Chosen for no particularly good reason.

            if( ( IITER / 20 )*20 == IITER && BSCALE*ABS1(T( ILAST, ILAST )) > SAFMIN ) {
               ESHIFT = ESHIFT + ( ASCALE*H( ILAST, ILAST ) )/( BSCALE*T( ILAST, ILAST ) );
            } else {
               ESHIFT = ESHIFT + ( ASCALE*H( ILAST, ILAST-1 ) )/( BSCALE*T( ILAST-1, ILAST-1 ) );
            }
            SHIFT = ESHIFT;
         }

         // Now check for two consecutive small subdiagonals.

         for (J = ILAST - 1; J >= IFIRST + 1; J--) { // 80
            ISTART = J;
            CTEMP = ASCALE*H( J, J ) - SHIFT*( BSCALE*T( J, J ) );
            TEMP = ABS1( CTEMP );
            TEMP2 = ASCALE*ABS1( H( J+1, J ) );
            TEMPR = max( TEMP, TEMP2 );
            if ( TEMPR < ONE && TEMPR != ZERO ) {
               TEMP = TEMP / TEMPR;
               TEMP2 = TEMP2 / TEMPR;
            }
            if( ABS1( H( J, J-1 ) )*TEMP2 <= TEMP*ATOL ) GO TO 90;
         } // 80

         ISTART = IFIRST;
         CTEMP = ASCALE*H( IFIRST, IFIRST ) - SHIFT*( BSCALE*T( IFIRST, IFIRST ) );
         } // 90

         // Do an implicit-shift QZ sweep.

         // Initial Q

         CTEMP2 = ASCALE*H( ISTART+1, ISTART );
         zlartg(CTEMP, CTEMP2, C, S, CTEMP3 );

         // Sweep

         for (J = ISTART; J <= ILAST - 1; J++) { // 150
            if ( J > ISTART ) {
               CTEMP = H( J, J-1 );
               zlartg(CTEMP, H( J+1, J-1 ), C, S, H( J, J-1 ) );
               H[J+1, J-1] = CZERO;
            }

            for (JC = J; JC <= ILASTM; JC++) { // 100
               CTEMP = C*H( J, JC ) + S*H( J+1, JC );
               H[J+1, JC] = -DCONJG( S )*H( J, JC ) + C*H( J+1, JC );
               H[J, JC] = CTEMP;
               CTEMP2 = C*T( J, JC ) + S*T( J+1, JC );
               T[J+1, JC] = -DCONJG( S )*T( J, JC ) + C*T( J+1, JC );
               T[J, JC] = CTEMP2;
            } // 100
            if ( ILQ ) {
               for (JR = 1; JR <= N; JR++) { // 110
                  CTEMP = C*Q( JR, J ) + DCONJG( S )*Q( JR, J+1 );
                  Q[JR, J+1] = -S*Q( JR, J ) + C*Q( JR, J+1 );
                  Q[JR, J] = CTEMP;
               } // 110
            }

            CTEMP = T( J+1, J+1 );
            zlartg(CTEMP, T( J+1, J ), C, S, T( J+1, J+1 ) );
            T[J+1, J] = CZERO;

            for (JR = IFRSTM; JR <= min( J+2, ILAST ); JR++) { // 120
               CTEMP = C*H( JR, J+1 ) + S*H( JR, J );
               H[JR, J] = -DCONJG( S )*H( JR, J+1 ) + C*H( JR, J );
               H[JR, J+1] = CTEMP;
            } // 120
            for (JR = IFRSTM; JR <= J; JR++) { // 130
               CTEMP = C*T( JR, J+1 ) + S*T( JR, J );
               T[JR, J] = -DCONJG( S )*T( JR, J+1 ) + C*T( JR, J );
               T[JR, J+1] = CTEMP;
            } // 130
            if ( ILZ ) {
               for (JR = 1; JR <= N; JR++) { // 140
                  CTEMP = C*Z( JR, J+1 ) + S*Z( JR, J );
                  Z[JR, J] = -DCONJG( S )*Z( JR, J+1 ) + C*Z( JR, J );
                  Z[JR, J+1] = CTEMP;
               } // 140
            }
         } // 150

         } // 160

      } // 170

      // Drop-through = non-convergence

      } // 180
      INFO = ILAST;
      GO TO 210;

      // Successful completion of all QZ steps

      } // 190

      // Set Eigenvalues 1:ILO-1

      for (J = 1; J <= ILO - 1; J++) { // 200
         ABSB = ( T( J, J ) ).abs();
         if ( ABSB > SAFMIN ) {
            SIGNBC = DCONJG( T( J, J ) / ABSB );
            T[J, J] = ABSB;
            if ( ILSCHR ) {
               zscal(J-1, SIGNBC, T( 1, J ), 1 );
               zscal(J, SIGNBC, H( 1, J ), 1 );
            } else {
               zscal(1, SIGNBC, H( J, J ), 1 );
            }
            if (ILZ) zscal( N, SIGNBC, Z( 1, J ), 1 );
         } else {
            T[J, J] = CZERO;
         }
         ALPHA[J] = H( J, J );
         BETA[J] = T( J, J );
      } // 200

      // Normal Termination

      INFO = 0;

      // Exit (other than argument error) -- return optimal workspace size

      } // 210
      WORK[1] = DCMPLX( N );
      return;
      }