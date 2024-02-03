
*>
*  =====================================================================
      SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, LDA, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             A( LDA, * ), W( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..

      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, INDE, INDTAU, INDWK2, INDWRK, ISCALE, LIOPT, LIWMIN, LLWORK, LLWRK2, LOPT, LWMIN;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, DLANSY;
      // EXTERNAL LSAME, DLAMCH, DLANSY, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACPY, DLASCL, DORMTR, DSCAL, DSTEDC, DSTERF, DSYTRD, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )

      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF

      IF( INFO.EQ.0 ) THEN
         IF( N.LE.1 ) THEN
            LIWMIN = 1
            LWMIN = 1
            LOPT = LWMIN
            LIOPT = LIWMIN
         } else {
            IF( WANTZ ) THEN
               LIWMIN = 3 + 5*N
               LWMIN = 1 + 6*N + 2*N**2
            } else {
               LIWMIN = 1
               LWMIN = 2*N + 1
            END IF
            LOPT = MAX( LWMIN, 2*N + N*ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 ) )
            LIOPT = LIWMIN
         END IF
         WORK( 1 ) = LOPT
         IWORK( 1 ) = LIOPT

         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -8
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -10
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYEVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      IF( N.EQ.1 ) THEN
         W( 1 ) = A( 1, 1 )
         IF( WANTZ ) A( 1, 1 ) = ONE
         RETURN
      END IF

      // Get machine constants.

      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )

      // Scale matrix to allowable range, if necessary.

      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) CALL DLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )

      // Call DSYTRD to reduce symmetric matrix to tridiagonal form.

      INDE = 1
      INDTAU = INDE + N
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      INDWK2 = INDWRK + N*N
      LLWRK2 = LWORK - INDWK2 + 1

      CALL DSYTRD( UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO )

      // For eigenvalues only, call DSTERF.  For eigenvectors, first call
      // DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
     t // ridiagonal matrix, then call DORMTR to multiply it by the
      // Householder transformations stored in A.

      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, WORK( INDE ), INFO )
      } else {
         CALL DSTEDC( 'I', N, W, WORK( INDE ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IWORK, LIWORK, INFO )          CALL DORMTR( 'L', UPLO, 'N', N, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IINFO )
         CALL DLACPY( 'A', N, N, WORK( INDWRK ), N, A, LDA )
      END IF

      // If matrix was scaled, then rescale eigenvalues appropriately.

      IF( ISCALE.EQ.1 ) CALL DSCAL( N, ONE / SIGMA, W, 1 )

      WORK( 1 ) = LOPT
      IWORK( 1 ) = LIOPT

      RETURN

      // End of DSYEVD

      }
