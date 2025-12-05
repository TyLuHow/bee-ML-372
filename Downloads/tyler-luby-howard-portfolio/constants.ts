import { Project, Experience, Resource, BlogPost, Paper } from './types';

// Images are placeholders as requested
const PLACEHOLDER_IMG = (w: number, h: number) => `https://picsum.photos/${w}/${h}`;

export const PROJECTS: Project[] = [
  {
    id: 'stormwater',
    title: 'Stormwater Watch',
    category: 'Web App',
    shortDescription: 'Real-time environmental monitoring platform tracking industrial violations.',
    fullDescription: 'A comprehensive dashboard designed for environmental compliance, tracking over 1.2M measurements across 369 facilities. The system aggregates real-time sensor data to detect and flag violations automatically.',
    metrics: ['1.2M+ monitored measurements', '369 facilities tracked', '899 parameters analyzed'],
    techStack: ['Next.js', 'Tailwind CSS', 'D3.js', 'PostgreSQL'],
    image: PLACEHOLDER_IMG(800, 600),
    link: 'https://stormwater-watch.vercel.app/dashboard',
  },
  {
    id: 'onroute',
    title: 'OnRoute Outdoors',
    category: 'Web App',
    shortDescription: 'Discovers outdoor activities along travel routes to break up long drives.',
    fullDescription: 'Solves the "missed adventure" problem by surfacing outdoor opportunities along a driving route. Features a custom routing algorithm that identifies points of interest within a configurable detour time.',
    metrics: ['Google Maps API Integration', 'Real-time routing', '100+ activity types'],
    techStack: ['Next.js', 'React Query', 'Google Maps API'],
    image: '/images/projects/onroute-outdoors/hero.png',
    link: 'https://onroute-outdoors.vercel.app/',
  },
  {
    id: 'hypertrack',
    title: 'HyperTrack Pro',
    category: 'Web App',
    shortDescription: 'Evidence-based workout tracking with progressive overload algorithms.',
    fullDescription: 'A mobile-first progressive web app focused on strength training data. It moves beyond simple logging to suggest weight increases based on previous performance and recovery metrics.',
    metrics: ['95+ Lighthouse Performance', '46+ Exercise Library', 'Zero-latency offline mode'],
    techStack: ['Next.js', 'TypeScript', 'Supabase', 'PWA'],
    image: PLACEHOLDER_IMG(800, 602),
    link: 'https://hypertrack-wheat.vercel.app/',
  },
  {
    id: 'beeml',
    title: 'Bee ML',
    category: 'ML/AI',
    shortDescription: 'Predictive modeling for bee toxicity in agricultural environments.',
    fullDescription: 'Leveraging machine learning to predict the toxicity levels of various agrochemicals on bee populations. This tool aids in making environmentally conscious decisions in farming.',
    metrics: ['92% Prediction Accuracy', 'Environmental Impact Focus'],
    techStack: ['Python', 'Scikit-Learn', 'FastAPI', 'React'],
    image: PLACEHOLDER_IMG(800, 603),
    link: 'https://bee-ml-372.vercel.app/',
  },
  {
    id: 'hatch',
    title: 'Hatch AI Integration',
    category: 'Enterprise',
    shortDescription: 'MCP servers and RAG knowledge systems for enterprise data.',
    fullDescription: 'Architected and deployed Multi-Context Protocol (MCP) servers to unify disparate data sources. Built a RAG system to allow natural language querying of internal documentation.',
    metrics: ['2.5M+ monthly records processed', 'Reduced query time by 60%'],
    techStack: ['LangChain', 'Vector DB', 'Python', 'Docker'],
    image: PLACEHOLDER_IMG(800, 604),
  },
  {
    id: 'rivian',
    title: 'Rivian Supply Chain ML',
    category: 'Enterprise',
    shortDescription: 'ML dashboards and XGBoost models for cost avoidance.',
    fullDescription: 'Developed predictive models to identify supply chain bottlenecks before they occurred. Created interactive Databricks dashboards for executive reporting.',
    metrics: ['$3M+ cost avoidance', 'Improved forecast accuracy by 15%'],
    techStack: ['Databricks', 'XGBoost', 'Hex', 'Python'],
    image: PLACEHOLDER_IMG(800, 605),
  },
  {
    id: 'tesla',
    title: 'Tesla Operations Automation',
    category: 'Enterprise',
    shortDescription: 'PowerApps automation for daily operations and reporting.',
    fullDescription: 'Streamlined manual reporting processes using Power Platform. Automated data ingestion for the daily operations review, saving significant analyst hours.',
    metrics: ['$35.7K quarterly savings', '400+ daily ops automated'],
    techStack: ['PowerApps', 'Power BI', 'SQL', 'DAX'],
    image: PLACEHOLDER_IMG(800, 606),
  },
  {
    id: 'rfnr',
    title: 'RFNR Accounting Hub',
    category: 'Enterprise',
    shortDescription: 'PostgreSQL hub integrating 10+ financial systems.',
    fullDescription: 'Designed a centralized data hub to synchronize financial data from various legacy systems. Built a GraphQL API layer for frontend consumption.',
    metrics: ['10+ Systems Integrated', 'Real-time reconciliation'],
    techStack: ['PostgreSQL', 'GraphQL', 'Node.js'],
    image: PLACEHOLDER_IMG(800, 607),
  }
];

export const EXPERIENCE: Experience[] = [
  {
    id: 'hatch-exp',
    company: 'Hatch',
    title: 'AI Engineer Intern',
    period: '2023 - Present',
    location: 'Remote',
    achievements: [
      'Developed MCP servers to standardize interface between LLMs and internal tools.',
      'Built an AI-powered ETL pipeline processing 2.5M+ monthly records.',
      'Implemented RAG systems for efficient internal knowledge retrieval.'
    ],
    techUsed: ['Python', 'LangChain', 'Pinecone', 'AWS']
  },
  {
    id: 'rivian-exp',
    company: 'Rivian',
    title: 'Supply Chain Data Scientist',
    period: '2022 - 2023',
    location: 'Palo Alto, CA',
    achievements: [
      'Deployed XGBoost models for demand forecasting, resulting in $1M-$3M cost avoidance.',
      'Created interactive Hex notebooks for cross-team data visualization.',
      'Optimized inventory allocation algorithms using Python and SQL.'
    ],
    techUsed: ['Databricks', 'Python', 'SQL', 'Tableau']
  },
  {
    id: 'tesla-exp',
    company: 'Tesla',
    title: 'Industrial Engineering Analyst',
    period: '2021 - 2022',
    location: 'Fremont, CA',
    achievements: [
      'Automated 400+ daily operations reports using Power Automate.',
      'Designed Power BI dashboards that identified $35.7K in quarterly savings.',
      'Conducted time studies and line balancing for Model 3 assembly.'
    ],
    techUsed: ['Power BI', 'PowerApps', 'Excel VBA', 'SQL']
  }
];

export const RESOURCES: Resource[] = [
  { id: '1', title: 'LangChain Cheatsheet', type: 'Tool', description: 'Quick reference for common chains and agents.', link: '#' },
  { id: '2', title: 'RAG Implementation Guide', type: 'Tutorial', description: 'Step-by-step production RAG setup.', link: '#' },
  { id: '3', title: 'Vector Search Demo', type: 'Demo', description: 'Interactive visualization of embeddings.', link: '#' },
  { id: '4', title: 'Prompt Engineering 101', type: 'Tutorial', description: 'Techniques for consistent LLM outputs.', link: '#' },
];

export const PAPERS: Paper[] = [
  { 
    id: 'p1', 
    title: 'Optimization of Supply Chain Networks via ML', 
    publication: 'IISE Annual Conference', 
    abstract: 'A comparative study of traditional OR methods vs. modern ML approaches in stochastic supply environments.',
    link: '#' 
  },
  { 
    id: 'p2', 
    title: 'Quantifying Soft Savings in Automation', 
    publication: 'Cal Poly Graduate Review', 
    abstract: 'Frameworks for measuring the financial impact of employee time savings through low-code automation.',
    link: '#' 
  }
];

export const BLOG_POSTS: BlogPost[] = [
  {
    id: 'b1',
    title: 'Bridging the Gap: IE to Software Engineering',
    excerpt: 'How industrial engineering principles apply directly to scalable software architecture.',
    date: 'Oct 12, 2023',
    readTime: '5 min read',
    tags: ['Career', 'Engineering']
  },
  {
    id: 'b2',
    title: 'Building a RAG System from Scratch',
    excerpt: 'Lessons learned implementing Retrieval Augmented Generation for enterprise data.',
    date: 'Sep 28, 2023',
    readTime: '8 min read',
    tags: ['AI', 'Tutorial']
  },
  {
    id: 'b3',
    title: 'The Modern Data Stack for Manufacturing',
    excerpt: 'Why factories need more than just Excel and legacy ERPs.',
    date: 'Aug 15, 2023',
    readTime: '6 min read',
    tags: ['Data', 'Industry 4.0']
  }
];
