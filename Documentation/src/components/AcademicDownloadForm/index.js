import React, {useEffect, useRef, useState} from 'react';
import styles from './styles.module.css';

const DOWNLOAD_PROCESSOR =
  'https://www.controlengineering.co.nz/assets/php/process_download_form.php';

// This is the public client identifier rendered by the original form.
// The corresponding verification secret remains only on the download server.
const PUBLIC_RECAPTCHA_SITE_KEY = '6Lc9bMQkAAAAAFX77g1d9g2iG31HPc-B6TzsAFtb';
const RECAPTCHA_SCRIPT_ID = 'opti-academic-recaptcha';

export default function AcademicDownloadForm() {
  const captchaContainer = useRef(null);
  const [captchaState, setCaptchaState] = useState('loading');

  useEffect(() => {
    let cancelled = false;
    let widgetId;

    const renderCaptcha = () => {
      if (
        cancelled ||
        !captchaContainer.current ||
        captchaContainer.current.dataset.rendered === 'true' ||
        !window.grecaptcha?.render
      ) {
        return;
      }

      window.grecaptcha.ready(() => {
        if (cancelled || !captchaContainer.current) {
          return;
        }
        widgetId = window.grecaptcha.render(captchaContainer.current, {
          sitekey: PUBLIC_RECAPTCHA_SITE_KEY,
          theme: document.documentElement.dataset.theme === 'dark' ? 'dark' : 'light',
        });
        captchaContainer.current.dataset.rendered = 'true';
        setCaptchaState('ready');
      });
    };

    let script = document.getElementById(RECAPTCHA_SCRIPT_ID);
    if (window.grecaptcha?.render) {
      renderCaptcha();
    } else if (script) {
      script.addEventListener('load', renderCaptcha);
    } else {
      script = document.createElement('script');
      script.id = RECAPTCHA_SCRIPT_ID;
      script.src = 'https://www.google.com/recaptcha/api.js?render=explicit';
      script.async = true;
      script.defer = true;
      script.addEventListener('load', renderCaptcha);
      script.addEventListener('error', () => {
        if (!cancelled) {
          setCaptchaState('error');
        }
      });
      document.head.appendChild(script);
    }

    return () => {
      cancelled = true;
      script?.removeEventListener('load', renderCaptcha);
      if (widgetId !== undefined && window.grecaptcha?.reset) {
        window.grecaptcha.reset(widgetId);
      }
    };
  }, []);

  return (
    <section className={styles.panel} aria-labelledby="academic-request-heading">
      <div className={styles.header}>
        <p className={styles.eyebrow}>Academic use only</p>
        <h2 id="academic-request-heading">Request the academic extension</h2>
        <p>
          All fields are required. Your request is submitted securely to the
          existing Control Engineering download service.
        </p>
      </div>

      <form action={DOWNLOAD_PROCESSOR} method="post" className={styles.form}>
        <input type="hidden" name="Software" value="optiAcademic" />

        <div className={styles.fieldGrid}>
          <div className={styles.field}>
            <label htmlFor="academic-full-name">Name</label>
            <input
              id="academic-full-name"
              name="FullName"
              type="text"
              autoComplete="name"
              required
            />
          </div>

          <div className={styles.field}>
            <label htmlFor="academic-email">Email</label>
            <input
              id="academic-email"
              name="EmailAddr"
              type="email"
              autoComplete="email"
              required
            />
          </div>

          <div className={`${styles.field} ${styles.fullWidth}`}>
            <label htmlFor="academic-institution">University / Institution</label>
            <input
              id="academic-institution"
              name="Company"
              type="text"
              autoComplete="organization"
              required
            />
          </div>

          <div className={`${styles.field} ${styles.fullWidth}`}>
            <label htmlFor="academic-intended-use">Intended use</label>
            <textarea
              id="academic-intended-use"
              name="IntendedUse"
              rows="4"
              required
            />
          </div>
        </div>

        <label className={styles.certification}>
          <input type="checkbox" name="DownloadR" value="no" required />
          <span>
            I certify that I will use this software only as a member of a
            noncommercial and academic institute.
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

        <button type="submit" className="button button--primary button--lg">
          Download
        </button>
      </form>
    </section>
  );
}
