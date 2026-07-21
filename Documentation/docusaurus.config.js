const lightCodeTheme = require('prism-react-renderer').themes.github;
const darkCodeTheme = require('prism-react-renderer').themes.dracula;
const migrationManifest = require('./migration/manifest.json');

const normalizeRoute = (route) => {
  if (route === '/') {
    return route;
  }
  return `/${route.replace(/^\/|\/$/g, '')}/`;
};

const legacyPageByRoute = new Map(
  migrationManifest.pages.map((page) => [
    normalizeRoute(page.route),
    page.source.replace('.', '/'),
  ]),
);

async function createConfig() {
  const remarkMath = (await import('remark-math')).default;
  const rehypeKatex = (await import('rehype-katex')).default;

  return {
    title: 'OPTI Toolbox',
    tagline:
      'A free MATLAB toolbox for constructing and solving optimization problems',
    url: process.env.SITE_URL || 'https://jonathancurrie.github.io',
    baseUrl: process.env.BASE_URL || '/OPTI/',
    organizationName: 'jonathancurrie',
    projectName: 'OPTI',
    trailingSlash: true,
    onBrokenLinks: 'throw',
    onBrokenAnchors: 'throw',
    onDuplicateRoutes: 'throw',

    markdown: {
      format: 'detect',
    },

    presets: [
      [
        'classic',
        {
          docs: {
            routeBasePath: '/',
            sidebarPath: require.resolve('./sidebars.js'),
            breadcrumbs: true,
            remarkPlugins: [remarkMath],
            rehypePlugins: [rehypeKatex],
          },
          blog: false,
          pages: false,
          theme: {
            customCss: require.resolve('./src/css/custom.css'),
          },
          gtag: {
            trackingID: 'G-9GBVJ653GH',
          },
          sitemap: {
            changefreq: 'monthly',
            priority: 0.5,
          },
        },
      ],
    ],

    plugins: [
      [
        '@docusaurus/plugin-client-redirects',
        {
          createRedirects(existingPath) {
            const legacyPage = legacyPageByRoute.get(normalizeRoute(existingPath));
            return legacyPage
              ? [
                  `/pmwiki.php/${legacyPage}`,
                  `/index.php/${legacyPage}`,
                ]
              : undefined;
          },
        },
      ],
      [
        require.resolve('@cmfcmf/docusaurus-search-local'),
        {
          indexDocs: true,
          indexBlog: false,
          indexPages: false,
          language: 'en',
          maxSearchResults: 10,
        },
      ],
    ],

    themeConfig: {
      metadata: [
        {
          name: 'description',
          content:
            'Documentation for OPTI Toolbox, a free MATLAB optimization toolbox.',
        },
      ],
      colorMode: {
        defaultMode: 'light',
        respectPrefersColorScheme: true,
      },
      navbar: {
        title: 'OPTI Toolbox',
        hideOnScroll: true,
        items: [
          {
            type: 'docSidebar',
            sidebarId: 'docsSidebar',
            position: 'left',
            label: 'Documentation',
          },
          {
            to: '/examples/',
            label: 'Examples',
            position: 'left',
          },
          {
            to: '/solvers/',
            label: 'Solvers',
            position: 'left',
          },
          {
            href: 'https://github.com/jonathancurrie/OPTI',
            label: 'GitHub',
            position: 'right',
          },
        ],
      },
      footer: {
        style: 'dark',
        links: [
          {
            title: 'Documentation',
            items: [
              {label: 'Getting started', to: '/getting-started/basics/'},
              {label: 'Examples', to: '/examples/'},
              {label: 'Solver reference', to: '/solvers/'},
            ],
          },
          {
            title: 'Project',
            items: [
              {label: 'Download', to: '/project/download/'},
              {label: 'License', to: '/project/license/'},
              {label: 'Citation', to: '/project/citation/'},
            ],
          },
          {
            title: 'Community',
            items: [
              {
                label: 'Question & Answer forum',
                href: 'https://groups.google.com/g/opti-toolbox-forum',
              },
              {
                label: 'OPTI on GitHub',
                href: 'https://github.com/jonathancurrie/OPTI',
              },
            ],
          },
        ],
        copyright: `Copyright © ${new Date().getFullYear()} <a href="https://www.controlengineering.co.nz/">Control Engineering</a>.`,
      },
      prism: {
        theme: lightCodeTheme,
        darkTheme: darkCodeTheme,
        additionalLanguages: ['matlab'],
      },
    },
  };
}

module.exports = createConfig;
