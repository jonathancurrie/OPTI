import React, {useEffect, useRef, useState} from 'react';
import styles from './styles.module.css';

const DOWNLOAD_URL =
  'https://www.dropbox.com/sh/7vtp9ifeuwt1h91/AACPsSpuJV8LBUhIImZu8pfCa?dl=0';
const PUBLIC_RECAPTCHA_SITE_KEY = '6Lfww10tAAAAAOBPqtuZWOmJlnwnD113ROgNIQhg';
const RECAPTCHA_ACTION = 'ACADEMIC_DOWNLOAD';
const RECAPTCHA_SCRIPT_ID = 'opti-academic-recaptcha-enterprise';
const RECAPTCHA_SCRIPT_URL =
  'https://www.google.com/recaptcha/enterprise.js?render=explicit';

export default function AcademicDownloadAccess() {
  const [academicUseConfirmed, setAcademicUseConfirmed] = useState(false);
  const [captchaState, setCaptchaState] = useState('loading');
  const captchaContainer = useRef(null);
  const downloadEnabled = academicUseConfirmed && captchaState === 'verified';

  useEffect(() => {
    let cancelled = false;
    let widgetId;

    const renderCaptcha = () => {
      const captcha = window.grecaptcha?.enterprise;
      if (
        cancelled ||
        !captchaContainer.current ||
        captchaContainer.current.dataset.rendered === 'true'
      ) {
        return;
      }

      if (!captcha?.render) {
        setCaptchaState('error');
        return;
      }

      captcha.ready(() => {
        if (cancelled || !captchaContainer.current) {
          return;
        }

        try {
          setCaptchaState('ready');
          widgetId = captcha.render(captchaContainer.current, {
            sitekey: PUBLIC_RECAPTCHA_SITE_KEY,
            action: RECAPTCHA_ACTION,
            theme:
              document.documentElement.dataset.theme === 'dark'
                ? 'dark'
                : 'light',
            callback: () => {
              if (!cancelled) {
                setCaptchaState('verified');
              }
            },
            'expired-callback': () => {
              if (!cancelled) {
                setCaptchaState('ready');
              }
            },
            'error-callback': () => {
              if (!cancelled) {
                setCaptchaState('error');
              }
            },
          });
          captchaContainer.current.dataset.rendered = 'true';
        } catch (error) {
          console.error('Unable to render reCAPTCHA Enterprise.', error);
          setCaptchaState('error');
        }
      });
    };

    const handleScriptError = () => {
      if (!cancelled) {
        setCaptchaState('error');
      }
    };

    let script = document.getElementById(RECAPTCHA_SCRIPT_ID);
    if (window.grecaptcha?.enterprise?.render) {
      renderCaptcha();
    } else if (script) {
      script.addEventListener('load', renderCaptcha);
      script.addEventListener('error', handleScriptError);
    } else {
      script = document.createElement('script');
      script.id = RECAPTCHA_SCRIPT_ID;
      script.src = RECAPTCHA_SCRIPT_URL;
      script.async = true;
      script.defer = true;
      script.addEventListener('load', renderCaptcha);
      script.addEventListener('error', handleScriptError);
      document.head.appendChild(script);
    }

    return () => {
      cancelled = true;
      script?.removeEventListener('load', renderCaptcha);
      script?.removeEventListener('error', handleScriptError);
      if (widgetId !== undefined && window.grecaptcha?.enterprise?.reset) {
        window.grecaptcha.enterprise.reset(widgetId);
      }
    };
  }, []);

  return (
    <section className={styles.panel} aria-labelledby="academic-download-heading">
      <div className={styles.header}>
        <p className={styles.eyebrow}>Academic use only</p>
        <h2 id="academic-download-heading">Access the academic extension</h2>
        <p>
          Confirm academic-only use and complete verification to enable the
          download link. This site does not request or store personal details.
        </p>
      </div>

      <div className={styles.accessControl}>
        <label className={styles.certification}>
          <input
            type="checkbox"
            checked={academicUseConfirmed}
            onChange={(event) => setAcademicUseConfirmed(event.target.checked)}
            aria-controls="academic-download-link"
          />
          <span>
            I confirm that the SCIP MEX file will only be used for academic
            purposes.
          </span>
        </label>

        <div className={styles.captchaRegion}>
          <div ref={captchaContainer} />
          {captchaState === 'loading' && (
            <p className={styles.status} role="status">
              Loading verification…
            </p>
          )}
          {captchaState === 'error' && (
            <p className={styles.error} role="alert">
              Verification could not be loaded. Check your connection and
              refresh the page.
            </p>
          )}
        </div>

        <div className={styles.linkRow}>
          {downloadEnabled ? (
            <a
              id="academic-download-link"
              className="button button--primary button--lg"
              href={DOWNLOAD_URL}
              target="_blank"
              rel="noopener noreferrer">
              SCIP MEX File Dropbox Link
            </a>
          ) : (
            <span
              id="academic-download-link"
              className={`button button--secondary button--lg ${styles.disabledLink}`}
              aria-disabled="true">
              SCIP MEX File Dropbox Link
            </span>
          )}
        </div>
      </div>
    </section>
  );
}
