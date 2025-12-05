export interface Project {
  id: string;
  title: string;
  category: 'Web App' | 'Enterprise' | 'ML/AI';
  shortDescription: string;
  fullDescription: string;
  metrics: string[];
  techStack: string[];
  image: string;
  link?: string;
  demoUrl?: string; // For iframe embedding
}

export interface Experience {
  id: string;
  company: string;
  title: string;
  period: string;
  location: string;
  achievements: string[];
  techUsed: string[];
}

export interface Resource {
  id: string;
  title: string;
  type: 'Tool' | 'Tutorial' | 'Demo';
  description: string;
  link: string;
}

export interface BlogPost {
  id: string;
  title: string;
  excerpt: string;
  date: string;
  readTime: string;
  tags: string[];
}

export interface Paper {
  id: string;
  title: string;
  publication: string;
  abstract: string;
  link: string;
}
