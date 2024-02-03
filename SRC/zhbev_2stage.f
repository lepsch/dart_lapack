      SUBROUTINE ZHBEV_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, RWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KD, LDAB, LDZ, N, LWORK;
      // ..
      // .. Array Arguments ..
      double             RWORK( * ), W( * );
      COMPLEX*16         AB( LDAB, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, WANTZ, LQUERY;
      int                IINFO, IMAX, INDE, INDWRK, INDRWK, ISCALE, LLWORK, LWMIN, LHTRD, LWTRD, IB, INDHOUS;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      double             DLAMCH, ZLANHB;
      // EXTERNAL LSAME, DLAMCH, ZLANHB, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSTERF, XERBLA, ZLASCL, ZSTEQR, ZHETRD_2STAGE, ZHETRD_HB2ST
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )

      INFO = 0
      if ( .NOT.( LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( KD.LT.0 ) {
         INFO = -4
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -6
      } else if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) {
         INFO = -9
      }

      if ( INFO.EQ.0 ) {
         if ( N.LE.1 ) {
            LWMIN = 1
            WORK( 1 ) = LWMIN
         } else {
            IB    = ILAENV2STAGE( 2, 'ZHETRD_HB2ST', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'ZHETRD_HB2ST', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'ZHETRD_HB2ST', JOBZ, N, KD, IB, -1 )
            LWMIN = LHTRD + LWTRD
            WORK( 1 )  = LWMIN
         }

         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) INFO = -11
      }

      if ( INFO.NE.0 ) {
         xerbla('ZHBEV_2STAGE ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( N.EQ.1 ) {
         if ( LOWER ) {
            W( 1 ) = DBLE( AB( 1, 1 ) )
         } else {
            W( 1 ) = DBLE( AB( KD+1, 1 ) )
         }
         IF( WANTZ ) Z( 1, 1 ) = ONE
         RETURN
      }

      // Get machine constants.

      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS    = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = SQRT( SMLNUM )
      RMAX   = SQRT( BIGNUM )

      // Scale matrix to allowable range, if necessary.

      ANRM = ZLANHB( 'M', UPLO, N, KD, AB, LDAB, RWORK )
      ISCALE = 0
      if ( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE.EQ.1 ) {
         if ( LOWER ) {
            zlascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         } else {
            zlascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         }
      }

      // Call ZHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form.

      INDE    = 1
      INDHOUS = 1
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1

      zhetrd_hb2st("N", JOBZ, UPLO, N, KD, AB, LDAB, W, RWORK( INDE ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEQR.

      if ( .NOT.WANTZ ) {
         dsterf(N, W, RWORK( INDE ), INFO );
      } else {
         INDRWK = INDE + N
         zsteqr(JOBZ, N, W, RWORK( INDE ), Z, LDZ, RWORK( INDRWK ), INFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if ( ISCALE.EQ.1 ) {
         if ( INFO.EQ.0 ) {
            IMAX = N
         } else {
            IMAX = INFO - 1
         }
         dscal(IMAX, ONE / SIGMA, W, 1 );
      }

      // Set WORK(1) to optimal workspace size.

      WORK( 1 ) = LWMIN

      RETURN

      // End of ZHBEV_2STAGE

      }
