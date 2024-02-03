      SUBROUTINE CTGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, JOB;
      int                INFO, LDA, LDB, LDVL, LDVR, LWORK, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      int                IWORK( * );
      REAL               DIF( * ), S( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      int                IDIFJB;
      const              ZERO = 0.0E+0, ONE = 1.0E+0, IDIFJB = 3 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, SOMCON, WANTBH, WANTDF, WANTS;
      int                I, IERR, IFST, ILST, K, KS, LWMIN, N1, N2;
      REAL               BIGNUM, COND, EPS, LNRM, RNRM, SCALE, SMLNUM
      COMPLEX            YHAX, YHBX
      // ..
      // .. Local Arrays ..
      COMPLEX            DUMMY( 1 ), DUMMY1( 1 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SCNRM2, SLAMCH, SLAPY2, SROUNDUP_LWORK
      COMPLEX            CDOTC
      // EXTERNAL LSAME, SCNRM2, SLAMCH, SLAPY2, SROUNDUP_LWORK, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CLACPY, CTGEXC, CTGSYL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, MAX
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTDF = LSAME( JOB, 'V' ) .OR. WANTBH

      SOMCON = LSAME( HOWMNY, 'S' )

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )

      if ( .NOT.WANTS .AND. .NOT.WANTDF ) {
         INFO = -1
      } else if ( .NOT.LSAME( HOWMNY, 'A' ) .AND. .NOT.SOMCON ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -8
      } else if ( WANTS .AND. LDVL.LT.N ) {
         INFO = -10
      } else if ( WANTS .AND. LDVR.LT.N ) {
         INFO = -12
      } else {

         // Set M to the number of eigenpairs for which condition numbers
         // are required, and test MM.

         if ( SOMCON ) {
            M = 0
            DO 10 K = 1, N
               IF( SELECT( K ) ) M = M + 1
   10       CONTINUE
         } else {
            M = N
         }

         if ( N.EQ.0 ) {
            LWMIN = 1
         } else if ( LSAME( JOB, 'V' ) .OR. LSAME( JOB, 'B' ) ) {
            LWMIN = 2*N*N
         } else {
            LWMIN = N
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)

         if ( MM.LT.M ) {
            INFO = -15
         } else if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -18
         }
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CTGSNA', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
      KS = 0
      DO 20 K = 1, N

         // Determine whether condition numbers are required for the k-th
         // eigenpair.

         if ( SOMCON ) {
            IF( .NOT.SELECT( K ) ) GO TO 20
         }

         KS = KS + 1

         if ( WANTS ) {

            // Compute the reciprocal condition number of the k-th
            // eigenvalue.

            RNRM = SCNRM2( N, VR( 1, KS ), 1 )
            LNRM = SCNRM2( N, VL( 1, KS ), 1 )
            CALL CGEMV( 'N', N, N, CMPLX( ONE, ZERO ), A, LDA, VR( 1, KS ), 1, CMPLX( ZERO, ZERO ), WORK, 1 )
            YHAX = CDOTC( N, WORK, 1, VL( 1, KS ), 1 )
            CALL CGEMV( 'N', N, N, CMPLX( ONE, ZERO ), B, LDB, VR( 1, KS ), 1, CMPLX( ZERO, ZERO ), WORK, 1 )
            YHBX = CDOTC( N, WORK, 1, VL( 1, KS ), 1 )
            COND = SLAPY2( ABS( YHAX ), ABS( YHBX ) )
            if ( COND.EQ.ZERO ) {
               S( KS ) = -ONE
            } else {
               S( KS ) = COND / ( RNRM*LNRM )
            }
         }

         if ( WANTDF ) {
            if ( N.EQ.1 ) {
               DIF( KS ) = SLAPY2( ABS( A( 1, 1 ) ), ABS( B( 1, 1 ) ) )
            } else {

               // Estimate the reciprocal condition number of the k-th
               // eigenvectors.

               // Copy the matrix (A, B) to the array WORK and move the
               // (k,k)th pair to the (1,1) position.

               CALL CLACPY( 'Full', N, N, A, LDA, WORK, N )
               CALL CLACPY( 'Full', N, N, B, LDB, WORK( N*N+1 ), N )
               IFST = K
               ILST = 1

               CALL CTGEXC( .FALSE., .FALSE., N, WORK, N, WORK( N*N+1 ), N, DUMMY, 1, DUMMY1, 1, IFST, ILST, IERR )

               if ( IERR.GT.0 ) {

                  // Ill-conditioned problem - swap rejected.

                  DIF( KS ) = ZERO
               } else {

                  // Reordering successful, solve generalized Sylvester
                  // equation for R and L,
                             // A22 * R - L * A11 = A12
                             // B22 * R - L * B11 = B12,
                  // and compute estimate of Difl[(A11,B11), (A22, B22)].

                  N1 = 1
                  N2 = N - N1
                  I = N*N + 1
                  CALL CTGSYL( 'N', IDIFJB, N2, N1, WORK( N*N1+N1+1 ), N, WORK, N, WORK( N1+1 ), N, WORK( N*N1+N1+I ), N, WORK( I ), N, WORK( N1+I ), N, SCALE, DIF( KS ), DUMMY, 1, IWORK, IERR )
               }
            }
         }

   20 CONTINUE
      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      RETURN

      // End of CTGSNA

      }
