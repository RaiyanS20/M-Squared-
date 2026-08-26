import { EventCategory } from '../domain/index.js';

/**
 * The starter dataset.
 *
 * Editorial rules for this project — read these before adding anything:
 *   1. FACTUAL ONLY. Every entry describes something that is documented and
 *      uncontroversial as to *whether* it happened.
 *   2. EVERY ENTRY CITES A SOURCE. The domain object refuses to be constructed
 *      without one, so this is enforced, not merely encouraged.
 *   3. NEUTRAL SUMMARIES. Describe what happened and why it mattered. Leave the
 *      judgement to the learner.
 *   4. An event is filed under the country whose territory or government it
 *      concerns. Events spanning many countries (a world war) appear once per
 *      country, described from that country's angle.
 *
 * Wikipedia is used as the citation because it is stable, free, and itself
 * carries references a learner can follow further. For a production history
 * product you would cite primary sources or an academic encyclopaedia.
 */

export interface SeedCountry {
  code: string;
  name: string;
  region: string;
}

export interface SeedEvent {
  countryCode: string;
  year: number;
  monthDay: string | null;
  title: string;
  summary: string;
  category: EventCategory;
  sourceUrl: string;
}

export const SEED_COUNTRIES: SeedCountry[] = [
  { code: 'AUS', name: 'Australia', region: 'Oceania' },
  { code: 'BRA', name: 'Brazil', region: 'South America' },
  { code: 'CHN', name: 'China', region: 'Asia' },
  { code: 'DEU', name: 'Germany', region: 'Europe' },
  { code: 'EGY', name: 'Egypt', region: 'Africa' },
  { code: 'FRA', name: 'France', region: 'Europe' },
  { code: 'GBR', name: 'United Kingdom', region: 'Europe' },
  { code: 'IND', name: 'India', region: 'Asia' },
  { code: 'JPN', name: 'Japan', region: 'Asia' },
  { code: 'RUS', name: 'Russia', region: 'Europe' },
  { code: 'USA', name: 'United States', region: 'North America' },
  { code: 'ZAF', name: 'South Africa', region: 'Africa' },
];

const wiki = (article: string): string => `https://en.wikipedia.org/wiki/${article}`;

export const SEED_EVENTS: SeedEvent[] = [
  // ---------------------------------------------------------------- United States
  {
    countryCode: 'USA', year: 1903, monthDay: '12-17',
    title: 'First powered aeroplane flight at Kitty Hawk',
    summary:
      'Orville and Wilbur Wright made the first sustained, controlled flight of a powered, heavier-than-air aircraft near Kitty Hawk, North Carolina. The first of four flights that day lasted 12 seconds and covered 120 feet.',
    category: EventCategory.Science, sourceUrl: wiki('Wright_Flyer'),
  },
  {
    countryCode: 'USA', year: 1906, monthDay: '04-18',
    title: 'San Francisco earthquake and fire',
    summary:
      'An earthquake estimated at magnitude 7.9 struck northern California, and the fires that followed burned for several days. Roughly 80 percent of San Francisco was destroyed and over half its population left homeless.',
    category: EventCategory.Disaster, sourceUrl: wiki('1906_San_Francisco_earthquake'),
  },
  {
    countryCode: 'USA', year: 1920, monthDay: '08-18',
    title: 'Nineteenth Amendment ratified, granting women the vote',
    summary:
      'Tennessee became the thirty-sixth state to ratify the Nineteenth Amendment, completing the three-quarters majority needed. It barred the states from denying the vote on the basis of sex, after a campaign lasting more than seventy years.',
    category: EventCategory.Society, sourceUrl: wiki('Nineteenth_Amendment_to_the_United_States_Constitution'),
  },
  {
    countryCode: 'USA', year: 1929, monthDay: '10-29',
    title: 'Wall Street Crash',
    summary:
      'Share prices on the New York Stock Exchange collapsed over several days, with the heaviest selling on 29 October, remembered as Black Tuesday. The crash marked the onset of the Great Depression in the United States.',
    category: EventCategory.Economy, sourceUrl: wiki('Wall_Street_Crash_of_1929'),
  },
  {
    countryCode: 'USA', year: 1941, monthDay: '12-07',
    title: 'Attack on Pearl Harbor',
    summary:
      'Japanese carrier aircraft attacked the US Pacific Fleet at Pearl Harbor, Hawaii, killing over 2,400 people. The United States declared war on Japan the following day, entering the Second World War.',
    category: EventCategory.Conflict, sourceUrl: wiki('Attack_on_Pearl_Harbor'),
  },
  {
    countryCode: 'USA', year: 1955, monthDay: '12-01',
    title: 'Rosa Parks arrested, beginning the Montgomery bus boycott',
    summary:
      'Rosa Parks was arrested in Montgomery, Alabama, for refusing to give up her bus seat to a white passenger. The 381-day boycott that followed brought Martin Luther King Jr. to national attention and ended with bus segregation ruled unconstitutional.',
    category: EventCategory.Society, sourceUrl: wiki('Montgomery_bus_boycott'),
  },
  {
    countryCode: 'USA', year: 1963, monthDay: '08-28',
    title: 'March on Washington for Jobs and Freedom',
    summary:
      'About a quarter of a million people gathered at the Lincoln Memorial to demand civil and economic rights for Black Americans, where Martin Luther King Jr. delivered the "I Have a Dream" speech. It helped build support for the Civil Rights Act of 1964.',
    category: EventCategory.Society, sourceUrl: wiki('March_on_Washington_for_Jobs_and_Freedom'),
  },
  {
    countryCode: 'USA', year: 1969, monthDay: '07-20',
    title: 'Apollo 11 lands the first humans on the Moon',
    summary:
      'Neil Armstrong and Buzz Aldrin landed the lunar module Eagle in the Sea of Tranquillity while Michael Collins orbited above. An estimated 600 million people watched the first moonwalk on television.',
    category: EventCategory.Science, sourceUrl: wiki('Apollo_11'),
  },
  {
    countryCode: 'USA', year: 2001, monthDay: '09-11',
    title: 'September 11 attacks',
    summary:
      'Four hijacked airliners were flown at targets in New York and Washington, destroying the World Trade Center towers and damaging the Pentagon; nearly 3,000 people were killed. The attacks reshaped US foreign and security policy for a generation.',
    category: EventCategory.Conflict, sourceUrl: wiki('September_11_attacks'),
  },
  {
    countryCode: 'USA', year: 2008, monthDay: '09-15',
    title: 'Lehman Brothers files for bankruptcy',
    summary:
      'The investment bank Lehman Brothers filed the largest bankruptcy in US history, with over 600 billion dollars in assets. Its collapse turned the subprime mortgage crisis into a global financial crisis.',
    category: EventCategory.Economy, sourceUrl: wiki('Bankruptcy_of_Lehman_Brothers'),
  },

  // -------------------------------------------------------------- United Kingdom
  {
    countryCode: 'GBR', year: 1901, monthDay: '01-22',
    title: 'Death of Queen Victoria ends the Victorian era',
    summary:
      'Queen Victoria died at Osborne House after a reign of 63 years, the longest of any British monarch to that point. Her son succeeded her as Edward VII, giving the Edwardian era its name.',
    category: EventCategory.Politics, sourceUrl: wiki('Queen_Victoria'),
  },
  {
    countryCode: 'GBR', year: 1912, monthDay: '04-15',
    title: 'RMS Titanic sinks on her maiden voyage',
    summary:
      'The British liner Titanic, which had sailed from Southampton, struck an iceberg in the North Atlantic and sank with the loss of about 1,500 lives. The disaster led directly to the International Convention for the Safety of Life at Sea.',
    category: EventCategory.Disaster, sourceUrl: wiki('Titanic'),
  },
  {
    countryCode: 'GBR', year: 1928, monthDay: '09-28',
    title: 'Alexander Fleming discovers penicillin',
    summary:
      'At St Mary\'s Hospital in London, Fleming noticed that a Penicillium mould contaminating a culture plate had killed the surrounding bacteria. Developed into a usable drug in Oxford during the 1940s, it became the first widely effective antibiotic.',
    category: EventCategory.Science, sourceUrl: wiki('Penicillin'),
  },
  {
    countryCode: 'GBR', year: 1936, monthDay: '12-11',
    title: 'Abdication of Edward VIII',
    summary:
      'Edward VIII gave up the throne rather than end his relationship with Wallis Simpson, an American who had been divorced twice. His brother became George VI, and the crisis tested the constitutional relationship between monarch, government and Dominions.',
    category: EventCategory.Politics, sourceUrl: wiki('Edward_VIII_abdication_crisis'),
  },
  {
    countryCode: 'GBR', year: 1940, monthDay: '07-10',
    title: 'Battle of Britain begins',
    summary:
      'The Luftwaffe opened a sustained air campaign against British shipping, airfields and cities, aiming to win the air superiority needed for an invasion. RAF Fighter Command held out, and the invasion was postponed indefinitely in September.',
    category: EventCategory.Conflict, sourceUrl: wiki('Battle_of_Britain'),
  },
  {
    countryCode: 'GBR', year: 1948, monthDay: '07-05',
    title: 'National Health Service founded',
    summary:
      'The NHS began operating under Health Minister Aneurin Bevan, offering healthcare free at the point of use and funded from general taxation. It remains one of the most significant institutions created by the post-war Labour government.',
    category: EventCategory.Society, sourceUrl: wiki('National_Health_Service'),
  },
  {
    countryCode: 'GBR', year: 1953, monthDay: '04-25',
    title: 'Structure of DNA published in Nature',
    summary:
      'James Watson and Francis Crick, working in Cambridge and drawing on X-ray diffraction data produced by Rosalind Franklin and Raymond Gosling, published the double-helix model of DNA. It became the foundation of modern molecular biology.',
    category: EventCategory.Science, sourceUrl: wiki('Molecular_Structure_of_Nucleic_Acids:_A_Structure_for_Deoxyribose_Nucleic_Acid'),
  },
  {
    countryCode: 'GBR', year: 1973, monthDay: '01-01',
    title: 'United Kingdom joins the European Economic Community',
    summary:
      'The UK became a member of the EEC alongside Ireland and Denmark, after two earlier applications had been vetoed by France. Membership was confirmed by referendum in 1975.',
    category: EventCategory.Politics, sourceUrl: wiki('1973_enlargement_of_the_European_Communities'),
  },
  {
    countryCode: 'GBR', year: 1979, monthDay: '05-04',
    title: 'Margaret Thatcher becomes the first woman Prime Minister',
    summary:
      'Thatcher took office after the Conservatives won the general election, and served until 1990. Her governments privatised state industries, curbed trade union power and reshaped the post-war economic consensus.',
    category: EventCategory.Politics, sourceUrl: wiki('Margaret_Thatcher'),
  },
  {
    countryCode: 'GBR', year: 2016, monthDay: '06-23',
    title: 'Referendum votes to leave the European Union',
    summary:
      'On a turnout of 72 percent, 51.9 percent voted to leave the EU. The result triggered the resignation of Prime Minister David Cameron and four years of withdrawal negotiations.',
    category: EventCategory.Politics, sourceUrl: wiki('2016_United_Kingdom_European_Union_membership_referendum'),
  },
  {
    countryCode: 'GBR', year: 2020, monthDay: '01-31',
    title: 'United Kingdom leaves the European Union',
    summary:
      'The UK formally ceased to be an EU member state after 47 years, entering a transition period that ran to the end of 2020. It was the first country to withdraw from the Union.',
    category: EventCategory.Politics, sourceUrl: wiki('Brexit'),
  },

  // ---------------------------------------------------------------------- France
  {
    countryCode: 'FRA', year: 1900, monthDay: '04-14',
    title: 'Exposition Universelle opens in Paris',
    summary:
      'The world\'s fair drew about 50 million visitors over seven months and introduced the Paris Métro, moving walkways and the Grand and Petit Palais. The 1900 Olympic Games were held alongside it.',
    category: EventCategory.Culture, sourceUrl: wiki('Exposition_Universelle_(1900)'),
  },
  {
    countryCode: 'FRA', year: 1903, monthDay: '07-01',
    title: 'First Tour de France',
    summary:
      'The newspaper L\'Auto organised a six-stage bicycle race around France to boost circulation; Maurice Garin won. It grew into the world\'s best-known cycling race.',
    category: EventCategory.Culture, sourceUrl: wiki('1903_Tour_de_France'),
  },
  {
    countryCode: 'FRA', year: 1905, monthDay: '12-09',
    title: 'Law on the Separation of the Churches and the State',
    summary:
      'The law ended state funding of religion and established French laïcité, guaranteeing freedom of worship while removing religion from public institutions. It remains a defining principle of the French Republic.',
    category: EventCategory.Politics, sourceUrl: wiki('1905_French_law_on_the_Separation_of_the_Churches_and_the_State'),
  },
  {
    countryCode: 'FRA', year: 1916, monthDay: '02-21',
    title: 'Battle of Verdun begins',
    summary:
      'German forces attacked the fortress city of Verdun in an attempt to bleed the French army white. The battle lasted almost ten months and caused an estimated 700,000 casualties on both sides.',
    category: EventCategory.Conflict, sourceUrl: wiki('Battle_of_Verdun'),
  },
  {
    countryCode: 'FRA', year: 1919, monthDay: '06-28',
    title: 'Treaty of Versailles signed',
    summary:
      'The treaty ending the First World War with Germany was signed in the Hall of Mirrors at Versailles, five years to the day after the assassination in Sarajevo. Its territorial and reparations terms shaped European politics between the wars.',
    category: EventCategory.Politics, sourceUrl: wiki('Treaty_of_Versailles'),
  },
  {
    countryCode: 'FRA', year: 1940, monthDay: '06-14',
    title: 'German forces enter Paris',
    summary:
      'Paris was declared an open city and occupied without a fight. An armistice followed on 22 June, dividing France into an occupied north and the Vichy-administered south.',
    category: EventCategory.Conflict, sourceUrl: wiki('Battle_of_France'),
  },
  {
    countryCode: 'FRA', year: 1944, monthDay: '08-25',
    title: 'Liberation of Paris',
    summary:
      'After an uprising by the French Resistance, the German garrison surrendered to Free French and American forces. Charles de Gaulle walked down the Champs-Élysées the following day.',
    category: EventCategory.Conflict, sourceUrl: wiki('Liberation_of_Paris'),
  },
  {
    countryCode: 'FRA', year: 1958, monthDay: '10-04',
    title: 'Constitution of the Fifth Republic promulgated',
    summary:
      'Approved by referendum the previous month, the new constitution greatly strengthened the presidency in response to the instability of the Fourth Republic and the Algerian crisis. De Gaulle became president in December.',
    category: EventCategory.Politics, sourceUrl: wiki('French_Fifth_Republic'),
  },
  {
    countryCode: 'FRA', year: 1968, monthDay: '05-13',
    title: 'May 1968 general strike',
    summary:
      'Student protests in Paris grew into a general strike involving roughly ten million workers, bringing the economy to a halt. The upheaval produced major wage agreements and lasting changes in French social attitudes.',
    category: EventCategory.Society, sourceUrl: wiki('May_1968_events_in_France'),
  },
  {
    countryCode: 'FRA', year: 1981, monthDay: '09-27',
    title: 'First TGV high-speed rail service opens',
    summary:
      'The Paris–Lyon TGV began commercial service at speeds up to 260 km/h, cutting the journey to about two hours. It launched the European high-speed rail network.',
    category: EventCategory.Science, sourceUrl: wiki('TGV'),
  },
  {
    countryCode: 'FRA', year: 2015, monthDay: '12-12',
    title: 'Paris Agreement on climate change adopted',
    summary:
      'Delegates from 196 parties agreed at COP21 to limit global warming to well below 2 °C above pre-industrial levels, with each country setting its own emissions targets. It replaced the Kyoto Protocol framework.',
    category: EventCategory.Politics, sourceUrl: wiki('Paris_Agreement'),
  },

  // --------------------------------------------------------------------- Germany
  {
    countryCode: 'DEU', year: 1918, monthDay: '11-09',
    title: 'German Revolution ends the monarchy',
    summary:
      'Amid naval mutiny and mass strikes, Kaiser Wilhelm II\'s abdication was announced and a republic proclaimed in Berlin. An armistice ending the First World War followed two days later.',
    category: EventCategory.Politics, sourceUrl: wiki('German_revolution_of_1918%E2%80%931919'),
  },
  {
    countryCode: 'DEU', year: 1919, monthDay: '08-11',
    title: 'Weimar Constitution comes into force',
    summary:
      'Germany\'s first democratic constitution introduced proportional representation, universal suffrage from age 20, and a directly elected president with emergency powers. Those emergency powers were later used to dismantle the republic.',
    category: EventCategory.Politics, sourceUrl: wiki('Weimar_Constitution'),
  },
  {
    countryCode: 'DEU', year: 1923, monthDay: '11-15',
    title: 'Hyperinflation ends with the Rentenmark',
    summary:
      'After prices rose to the point where a US dollar cost over four trillion marks, a new currency backed by a mortgage on German land stabilised the economy. The episode left a lasting German aversion to inflation.',
    category: EventCategory.Economy, sourceUrl: wiki('Hyperinflation_in_the_Weimar_Republic'),
  },
  {
    countryCode: 'DEU', year: 1933, monthDay: '01-30',
    title: 'Hitler appointed Chancellor',
    summary:
      'President Hindenburg appointed Adolf Hitler head of a coalition government. Within months the Reichstag Fire Decree and the Enabling Act had removed civil liberties and parliamentary checks, establishing a dictatorship.',
    category: EventCategory.Politics, sourceUrl: wiki('Machtergreifung'),
  },
  {
    countryCode: 'DEU', year: 1939, monthDay: '09-01',
    title: 'Invasion of Poland begins the Second World War in Europe',
    summary:
      'German forces crossed into Poland, and Britain and France declared war two days later. The Soviet Union invaded eastern Poland on 17 September under the Molotov–Ribbentrop Pact.',
    category: EventCategory.Conflict, sourceUrl: wiki('Invasion_of_Poland'),
  },
  {
    countryCode: 'DEU', year: 1945, monthDay: '05-08',
    title: 'Unconditional surrender ends the war in Europe',
    summary:
      'Germany\'s surrender took effect, ending nearly six years of war in Europe in which tens of millions died, including six million Jews murdered in the Holocaust. The country was divided into four occupation zones.',
    category: EventCategory.Conflict, sourceUrl: wiki('End_of_World_War_II_in_Europe'),
  },
  {
    countryCode: 'DEU', year: 1949, monthDay: '05-23',
    title: 'Federal Republic of Germany founded',
    summary:
      'The Basic Law came into force in the three western occupation zones, creating a federal parliamentary democracy with strong constitutional protections. The German Democratic Republic was established in the Soviet zone in October.',
    category: EventCategory.Politics, sourceUrl: wiki('Basic_Law_for_the_Federal_Republic_of_Germany'),
  },
  {
    countryCode: 'DEU', year: 1961, monthDay: '08-13',
    title: 'Construction of the Berlin Wall begins',
    summary:
      'East German authorities sealed the border around West Berlin overnight, ending an exodus that had taken some 2.7 million people west since 1949. The wall stood for 28 years.',
    category: EventCategory.Politics, sourceUrl: wiki('Berlin_Wall'),
  },
  {
    countryCode: 'DEU', year: 1989, monthDay: '11-09',
    title: 'Fall of the Berlin Wall',
    summary:
      'A confused announcement about relaxed travel rules brought crowds to the checkpoints, and guards opened the barriers. The night became the symbol of the end of communist rule in Eastern Europe.',
    category: EventCategory.Politics, sourceUrl: wiki('Fall_of_the_Berlin_Wall'),
  },
  {
    countryCode: 'DEU', year: 1990, monthDay: '10-03',
    title: 'German reunification',
    summary:
      'The German Democratic Republic acceded to the Federal Republic, restoring a single German state less than a year after the wall opened. The date is now the national holiday.',
    category: EventCategory.Politics, sourceUrl: wiki('German_reunification'),
  },
  {
    countryCode: 'DEU', year: 2002, monthDay: '01-01',
    title: 'Euro notes and coins replace the Deutsche Mark',
    summary:
      'Germany exchanged one of the world\'s most trusted currencies for the euro, alongside eleven other countries. The changeover was the largest cash conversion ever attempted.',
    category: EventCategory.Economy, sourceUrl: wiki('Euro'),
  },

  // ---------------------------------------------------------------------- Russia
  {
    countryCode: 'RUS', year: 1905, monthDay: '01-22',
    title: 'Bloody Sunday in St Petersburg',
    summary:
      'Troops fired on a peaceful procession of workers carrying a petition to the Winter Palace, killing scores of demonstrators. The killings destroyed popular faith in the Tsar and set off the 1905 Revolution.',
    category: EventCategory.Society, sourceUrl: wiki('Bloody_Sunday_(1905)'),
  },
  {
    countryCode: 'RUS', year: 1917, monthDay: '03-08',
    title: 'February Revolution topples the Tsar',
    summary:
      'Bread protests and strikes in Petrograd, dated February in the old Russian calendar, spread until the garrison mutinied and Nicholas II abdicated. A Provisional Government took power alongside the workers\' soviets.',
    category: EventCategory.Politics, sourceUrl: wiki('February_Revolution'),
  },
  {
    countryCode: 'RUS', year: 1917, monthDay: '11-07',
    title: 'October Revolution brings the Bolsheviks to power',
    summary:
      'Bolshevik forces seized government buildings in Petrograd and deposed the Provisional Government. The seizure of power led to civil war and to the creation of the world\'s first communist state.',
    category: EventCategory.Politics, sourceUrl: wiki('October_Revolution'),
  },
  {
    countryCode: 'RUS', year: 1922, monthDay: '12-30',
    title: 'Soviet Union formally established',
    summary:
      'A treaty joined the Russian, Ukrainian, Byelorussian and Transcaucasian republics into the Union of Soviet Socialist Republics. It would last 69 years.',
    category: EventCategory.Politics, sourceUrl: wiki('Soviet_Union'),
  },
  {
    countryCode: 'RUS', year: 1941, monthDay: '06-22',
    title: 'Operation Barbarossa: Germany invades the Soviet Union',
    summary:
      'Some three million Axis troops attacked along a front of more than 2,000 kilometres, breaking the 1939 non-aggression pact. The Eastern Front became the deadliest theatre of the war.',
    category: EventCategory.Conflict, sourceUrl: wiki('Operation_Barbarossa'),
  },
  {
    countryCode: 'RUS', year: 1943, monthDay: '02-02',
    title: 'Battle of Stalingrad ends',
    summary:
      'The surrender of the encircled German Sixth Army ended five months of fighting that cost well over a million casualties in total. It marked the turning point of the war on the Eastern Front.',
    category: EventCategory.Conflict, sourceUrl: wiki('Battle_of_Stalingrad'),
  },
  {
    countryCode: 'RUS', year: 1957, monthDay: '10-04',
    title: 'Sputnik 1 becomes the first artificial satellite',
    summary:
      'The Soviet Union launched an 84-kilogram sphere into low Earth orbit, where its radio beeps could be picked up worldwide. The launch began the space age and the space race.',
    category: EventCategory.Science, sourceUrl: wiki('Sputnik_1'),
  },
  {
    countryCode: 'RUS', year: 1961, monthDay: '04-12',
    title: 'Yuri Gagarin becomes the first human in space',
    summary:
      'Gagarin completed a single orbit of the Earth aboard Vostok 1 in 108 minutes. He returned safely by parachute after ejecting from the capsule.',
    category: EventCategory.Science, sourceUrl: wiki('Vostok_1'),
  },
  {
    countryCode: 'RUS', year: 1991, monthDay: '12-26',
    title: 'Dissolution of the Soviet Union',
    summary:
      'The Supreme Soviet formally dissolved the USSR, a day after Mikhail Gorbachev resigned as president. Fifteen independent states emerged, with the Russian Federation as the legal successor.',
    category: EventCategory.Politics, sourceUrl: wiki('Dissolution_of_the_Soviet_Union'),
  },
  {
    countryCode: 'RUS', year: 1998, monthDay: '08-17',
    title: 'Russian financial crisis and rouble default',
    summary:
      'Russia devalued the rouble and defaulted on domestic debt after falling oil prices and fiscal strain. Output contracted sharply, though a weaker currency helped the economy recover within two years.',
    category: EventCategory.Economy, sourceUrl: wiki('1998_Russian_financial_crisis'),
  },

  // ----------------------------------------------------------------------- China
  {
    countryCode: 'CHN', year: 1900, monthDay: '06-20',
    title: 'Siege of the foreign legations during the Boxer Uprising',
    summary:
      'Boxer militia and imperial troops besieged the diplomatic quarter in Peking for 55 days. An eight-nation expeditionary force lifted the siege in August and imposed a heavy indemnity on the Qing court.',
    category: EventCategory.Conflict, sourceUrl: wiki('Boxer_Rebellion'),
  },
  {
    countryCode: 'CHN', year: 1911, monthDay: '10-10',
    title: 'Wuchang Uprising starts the Xinhai Revolution',
    summary:
      'A mutiny by army units in Wuchang triggered declarations of independence across the provinces. Within months the Qing dynasty, and with it two thousand years of imperial rule, had ended.',
    category: EventCategory.Politics, sourceUrl: wiki('Wuchang_Uprising'),
  },
  {
    countryCode: 'CHN', year: 1912, monthDay: '01-01',
    title: 'Republic of China proclaimed',
    summary:
      'Sun Yat-sen was inaugurated as provisional president in Nanjing. The last emperor, Puyi, abdicated on 12 February.',
    category: EventCategory.Politics, sourceUrl: wiki('Republic_of_China_(1912%E2%80%931949)'),
  },
  {
    countryCode: 'CHN', year: 1919, monthDay: '05-04',
    title: 'May Fourth Movement',
    summary:
      'Students in Beijing protested the Treaty of Versailles transferring German concessions in Shandong to Japan. The movement fuelled a broader push for science, vernacular literature and political reform.',
    category: EventCategory.Society, sourceUrl: wiki('May_Fourth_Movement'),
  },
  {
    countryCode: 'CHN', year: 1934, monthDay: '10-16',
    title: 'Long March begins',
    summary:
      'Communist forces broke out of Nationalist encirclement in Jiangxi and began a retreat of some 9,000 kilometres to Shaanxi. Perhaps one in ten survived, and Mao Zedong emerged as the party\'s leader.',
    category: EventCategory.Conflict, sourceUrl: wiki('Long_March'),
  },
  {
    countryCode: 'CHN', year: 1937, monthDay: '07-07',
    title: 'Marco Polo Bridge Incident begins full-scale war with Japan',
    summary:
      'A clash near Beijing escalated into the Second Sino-Japanese War, which lasted until 1945 and cost millions of Chinese lives. It merged into the wider Second World War after 1941.',
    category: EventCategory.Conflict, sourceUrl: wiki('Marco_Polo_Bridge_Incident'),
  },
  {
    countryCode: 'CHN', year: 1949, monthDay: '10-01',
    title: 'People\'s Republic of China proclaimed',
    summary:
      'Mao Zedong declared the new state from Tiananmen Gate after the Communist victory in the civil war. The Nationalist government withdrew to Taiwan.',
    category: EventCategory.Politics, sourceUrl: wiki('Proclamation_of_the_People%27s_Republic_of_China'),
  },
  {
    countryCode: 'CHN', year: 1966, monthDay: '05-16',
    title: 'Cultural Revolution launched',
    summary:
      'A party circular set off a decade of political campaigns in which schools closed, officials were purged and cultural heritage was destroyed. It was formally repudiated by the party in 1981.',
    category: EventCategory.Politics, sourceUrl: wiki('Cultural_Revolution'),
  },
  {
    countryCode: 'CHN', year: 1978, monthDay: '12-18',
    title: 'Reform and opening up begins',
    summary:
      'The Third Plenum of the Eleventh Central Committee shifted priority from class struggle to economic modernisation under Deng Xiaoping. The reforms that followed lifted hundreds of millions out of poverty.',
    category: EventCategory.Economy, sourceUrl: wiki('Chinese_economic_reform'),
  },
  {
    countryCode: 'CHN', year: 2001, monthDay: '12-11',
    title: 'China joins the World Trade Organization',
    summary:
      'Accession after fifteen years of negotiation bound China to lower tariffs and opened export markets. Chinese trade volumes multiplied over the following decade.',
    category: EventCategory.Economy, sourceUrl: wiki('China_and_the_World_Trade_Organization'),
  },
  {
    countryCode: 'CHN', year: 2008, monthDay: '08-08',
    title: 'Beijing Summer Olympics open',
    summary:
      'China hosted the Games for the first time, spending heavily on venues and infrastructure, and topped the gold medal table. The opening ceremony was watched by an estimated one billion people.',
    category: EventCategory.Culture, sourceUrl: wiki('2008_Summer_Olympics'),
  },

  // ----------------------------------------------------------------------- India
  {
    countryCode: 'IND', year: 1905, monthDay: '10-16',
    title: 'Partition of Bengal',
    summary:
      'The colonial government split Bengal into two provinces, a decision widely seen as designed to divide the nationalist movement. Mass protests and the Swadeshi boycott of British goods followed, and the partition was reversed in 1911.',
    category: EventCategory.Politics, sourceUrl: wiki('Partition_of_Bengal_(1905)'),
  },
  {
    countryCode: 'IND', year: 1919, monthDay: '04-13',
    title: 'Jallianwala Bagh massacre',
    summary:
      'Troops under Colonel Dyer fired without warning on a large gathering in an enclosed garden in Amritsar; official figures counted 379 dead, Indian estimates far more. The massacre turned moderate opinion decisively against British rule.',
    category: EventCategory.Conflict, sourceUrl: wiki('Jallianwala_Bagh_massacre'),
  },
  {
    countryCode: 'IND', year: 1930, monthDay: '03-12',
    title: 'Salt March begins',
    summary:
      'Gandhi walked 385 kilometres from Sabarmati to Dandi to make salt from seawater in defiance of the salt tax. The march made civil disobedience a mass movement and drew worldwide attention.',
    category: EventCategory.Society, sourceUrl: wiki('Salt_March'),
  },
  {
    countryCode: 'IND', year: 1942, monthDay: '08-08',
    title: 'Quit India Movement launched',
    summary:
      'The Congress demanded an immediate end to British rule; its leadership was arrested within hours and the movement suppressed. It nonetheless made continued British rule after the war politically untenable.',
    category: EventCategory.Politics, sourceUrl: wiki('Quit_India_Movement'),
  },
  {
    countryCode: 'IND', year: 1947, monthDay: '08-15',
    title: 'Independence and partition',
    summary:
      'India became independent as British India was partitioned into India and Pakistan. The accompanying migration displaced some 15 million people and communal violence killed hundreds of thousands.',
    category: EventCategory.Politics, sourceUrl: wiki('Partition_of_India'),
  },
  {
    countryCode: 'IND', year: 1950, monthDay: '01-26',
    title: 'Constitution of India comes into force',
    summary:
      'India became a republic under a constitution drafted by an assembly chaired in its drafting committee by B. R. Ambedkar. It remains the longest written constitution of any sovereign country.',
    category: EventCategory.Politics, sourceUrl: wiki('Constitution_of_India'),
  },
  {
    countryCode: 'IND', year: 1966, monthDay: '01-24',
    title: 'Indira Gandhi becomes Prime Minister',
    summary:
      'Gandhi took office as India\'s first woman prime minister, serving until 1977 and again from 1980. Her tenure included the 1971 war, the nuclear test of 1974 and the 1975–77 Emergency.',
    category: EventCategory.Politics, sourceUrl: wiki('Indira_Gandhi'),
  },
  {
    countryCode: 'IND', year: 1971, monthDay: '12-03',
    title: 'Indo-Pakistani War and the creation of Bangladesh',
    summary:
      'India intervened in the conflict in East Pakistan, and the war ended in thirteen days with the surrender of Pakistani forces in Dhaka. Bangladesh became independent.',
    category: EventCategory.Conflict, sourceUrl: wiki('Indo-Pakistani_War_of_1971'),
  },
  {
    countryCode: 'IND', year: 1991, monthDay: '07-24',
    title: 'Economic liberalisation begins',
    summary:
      'Facing a balance of payments crisis, Finance Minister Manmohan Singh presented a budget that dismantled industrial licensing, cut tariffs and opened the economy to foreign investment. Growth accelerated markedly in the decades that followed.',
    category: EventCategory.Economy, sourceUrl: wiki('Economic_liberalisation_in_India'),
  },
  {
    countryCode: 'IND', year: 2013, monthDay: '11-05',
    title: 'Mars Orbiter Mission launched',
    summary:
      'ISRO launched Mangalyaan, which entered Mars orbit in September 2014 at the first attempt on a budget of about 74 million dollars. India became the first Asian nation to reach Mars orbit.',
    category: EventCategory.Science, sourceUrl: wiki('Mars_Orbiter_Mission'),
  },

  // ----------------------------------------------------------------------- Japan
  {
    countryCode: 'JPN', year: 1904, monthDay: '02-08',
    title: 'Russo-Japanese War begins',
    summary:
      'Japan attacked the Russian fleet at Port Arthur over rival claims in Manchuria and Korea. Japan\'s victory in 1905 was the first modern defeat of a European power by an Asian state.',
    category: EventCategory.Conflict, sourceUrl: wiki('Russo-Japanese_War'),
  },
  {
    countryCode: 'JPN', year: 1923, monthDay: '09-01',
    title: 'Great Kantō earthquake',
    summary:
      'A magnitude 7.9 earthquake struck the Tokyo–Yokohama region at lunchtime, and the resulting firestorms killed over 100,000 people. Reconstruction reshaped Tokyo, and the disaster was followed by mob violence against Korean residents.',
    category: EventCategory.Disaster, sourceUrl: wiki('1923_Great_Kant%C5%8D_earthquake'),
  },
  {
    countryCode: 'JPN', year: 1945, monthDay: '08-06',
    title: 'Atomic bombing of Hiroshima',
    summary:
      'A US aircraft dropped an atomic bomb on Hiroshima, killing an estimated 70,000 people immediately and many more from injuries and radiation. Nagasaki was bombed on 9 August.',
    category: EventCategory.Conflict, sourceUrl: wiki('Atomic_bombings_of_Hiroshima_and_Nagasaki'),
  },
  {
    countryCode: 'JPN', year: 1945, monthDay: '08-15',
    title: 'Japan announces surrender',
    summary:
      'Emperor Hirohito announced acceptance of the Potsdam Declaration in a radio broadcast, the first time most citizens had heard his voice. Allied occupation under General MacArthur followed until 1952.',
    category: EventCategory.Conflict, sourceUrl: wiki('Surrender_of_Japan'),
  },
  {
    countryCode: 'JPN', year: 1947, monthDay: '05-03',
    title: 'Post-war Constitution takes effect',
    summary:
      'The new constitution made the emperor a symbol of the state, guaranteed civil liberties, and in Article 9 renounced war as a sovereign right. It has never been amended.',
    category: EventCategory.Politics, sourceUrl: wiki('Constitution_of_Japan'),
  },
  {
    countryCode: 'JPN', year: 1964, monthDay: '10-01',
    title: 'Tōkaidō Shinkansen opens',
    summary:
      'The world\'s first high-speed rail line linked Tokyo and Osaka at up to 210 km/h, nine days before the Tokyo Olympics opened. It became the model for high-speed rail everywhere.',
    category: EventCategory.Science, sourceUrl: wiki('Tokaido_Shinkansen'),
  },
  {
    countryCode: 'JPN', year: 1964, monthDay: '10-10',
    title: 'Tokyo Summer Olympics open',
    summary:
      'The first Olympics held in Asia showcased Japan\'s post-war recovery and were the first to be broadcast internationally by satellite. Colour television and computerised timing debuted at the Games.',
    category: EventCategory.Culture, sourceUrl: wiki('1964_Summer_Olympics'),
  },
  {
    countryCode: 'JPN', year: 1989, monthDay: '01-07',
    title: 'Death of Emperor Hirohito and the start of the Heisei era',
    summary:
      'Hirohito died after a 62-year reign spanning militarism, war, occupation and the economic miracle. His son Akihito succeeded him and the era name changed to Heisei.',
    category: EventCategory.Politics, sourceUrl: wiki('Hirohito'),
  },
  {
    countryCode: 'JPN', year: 1991, monthDay: null,
    title: 'Asset price bubble collapses, beginning the "lost decade"',
    summary:
      'Land and share prices fell sharply after the Bank of Japan tightened policy, ending the boom of the late 1980s. Growth and inflation stayed weak for many years afterwards.',
    category: EventCategory.Economy, sourceUrl: wiki('Japanese_asset_price_bubble'),
  },
  {
    countryCode: 'JPN', year: 2011, monthDay: '03-11',
    title: 'Tōhoku earthquake, tsunami and Fukushima accident',
    summary:
      'A magnitude 9.0 earthquake off north-east Japan generated a tsunami that killed around 18,000 people and caused meltdowns at the Fukushima Daiichi nuclear plant. It was the most powerful earthquake ever recorded in Japan.',
    category: EventCategory.Disaster, sourceUrl: wiki('2011_T%C5%8Dhoku_earthquake_and_tsunami'),
  },

  // ---------------------------------------------------------------------- Brazil
  {
    countryCode: 'BRA', year: 1917, monthDay: '10-26',
    title: 'Brazil declares war on Germany',
    summary:
      'After repeated sinkings of Brazilian merchant ships, Brazil entered the First World War, the only South American country to do so. Its contribution was mainly naval patrols and medical missions.',
    category: EventCategory.Conflict, sourceUrl: wiki('Brazil_during_World_War_I'),
  },
  {
    countryCode: 'BRA', year: 1922, monthDay: '02-13',
    title: 'Modern Art Week in São Paulo',
    summary:
      'A festival of exhibitions, concerts and readings at the Municipal Theatre launched Brazilian modernism. It set out to build a distinctly Brazilian art rather than imitate European models.',
    category: EventCategory.Culture, sourceUrl: wiki('Modern_Art_Week'),
  },
  {
    countryCode: 'BRA', year: 1930, monthDay: '10-24',
    title: 'Revolution of 1930 brings Getúlio Vargas to power',
    summary:
      'A military coup ended the Old Republic and installed Vargas, who governed for fifteen years, including the authoritarian Estado Novo. Labour laws and industrialisation date from this period.',
    category: EventCategory.Politics, sourceUrl: wiki('Brazilian_Revolution_of_1930'),
  },
  {
    countryCode: 'BRA', year: 1950, monthDay: '07-16',
    title: 'Maracanaço: Brazil loses the World Cup final at home',
    summary:
      'Uruguay beat Brazil 2–1 in front of a crowd officially recorded at about 174,000 at the Maracanã. The defeat became a lasting national reference point.',
    category: EventCategory.Culture, sourceUrl: wiki('Maracanazo'),
  },
  {
    countryCode: 'BRA', year: 1960, monthDay: '04-21',
    title: 'Brasília inaugurated as the new capital',
    summary:
      'Built in 41 months on the central plateau to designs by Lúcio Costa and Oscar Niemeyer, the city moved the capital inland from Rio de Janeiro. It is now a UNESCO World Heritage Site.',
    category: EventCategory.Politics, sourceUrl: wiki('Bras%C3%ADlia'),
  },
  {
    countryCode: 'BRA', year: 1964, monthDay: '03-31',
    title: 'Military coup begins 21 years of dictatorship',
    summary:
      'The armed forces deposed President João Goulart and governed until 1985, censoring the press and repressing opposition. The period also saw rapid, debt-financed economic growth.',
    category: EventCategory.Politics, sourceUrl: wiki('1964_Brazilian_coup_d%27%C3%A9tat'),
  },
  {
    countryCode: 'BRA', year: 1985, monthDay: '03-15',
    title: 'Return to civilian government',
    summary:
      'José Sarney was sworn in after the death of president-elect Tancredo Neves, ending military rule. Direct presidential elections returned in 1989.',
    category: EventCategory.Politics, sourceUrl: wiki('New_Republic_(Brazil)'),
  },
  {
    countryCode: 'BRA', year: 1988, monthDay: '10-05',
    title: 'New Constitution promulgated',
    summary:
      'Known as the "citizen constitution", it guaranteed broad social rights, strengthened the judiciary and decentralised power to states and municipalities. It remains in force.',
    category: EventCategory.Politics, sourceUrl: wiki('Constitution_of_Brazil'),
  },
  {
    countryCode: 'BRA', year: 1994, monthDay: '07-01',
    title: 'Plano Real ends hyperinflation',
    summary:
      'A new currency, the real, was introduced alongside fiscal and monetary reforms after years of inflation running above 2,000 percent annually. Inflation fell to single digits within a year.',
    category: EventCategory.Economy, sourceUrl: wiki('Plano_Real'),
  },
  {
    countryCode: 'BRA', year: 2016, monthDay: '08-05',
    title: 'Rio de Janeiro Summer Olympics open',
    summary:
      'Rio became the first South American city to host the Games. Preparations took place during a severe recession and a political crisis that led to the president\'s impeachment that month.',
    category: EventCategory.Culture, sourceUrl: wiki('2016_Summer_Olympics'),
  },

  // ---------------------------------------------------------------- South Africa
  {
    countryCode: 'ZAF', year: 1910, monthDay: '05-31',
    title: 'Union of South Africa established',
    summary:
      'Four British colonies were merged into a self-governing dominion. The franchise was restricted almost entirely to white voters, setting the pattern for the century that followed.',
    category: EventCategory.Politics, sourceUrl: wiki('Union_of_South_Africa'),
  },
  {
    countryCode: 'ZAF', year: 1912, monthDay: '01-08',
    title: 'Founding of the African National Congress',
    summary:
      'Delegates met in Bloemfontein to form the South African Native National Congress, renamed the ANC in 1923, to campaign for African political rights. It became the leading anti-apartheid organisation and later the governing party.',
    category: EventCategory.Politics, sourceUrl: wiki('African_National_Congress'),
  },
  {
    countryCode: 'ZAF', year: 1948, monthDay: '05-26',
    title: 'National Party wins election and begins apartheid',
    summary:
      'The National Party took office on a platform of apartheid, and over the following years passed laws classifying the population by race and segregating residence, marriage, work and education.',
    category: EventCategory.Politics, sourceUrl: wiki('Apartheid'),
  },
  {
    countryCode: 'ZAF', year: 1960, monthDay: '03-21',
    title: 'Sharpeville massacre',
    summary:
      'Police fired on a crowd protesting the pass laws, killing 69 people. The government declared a state of emergency and banned the ANC and PAC, pushing both towards armed struggle.',
    category: EventCategory.Conflict, sourceUrl: wiki('Sharpeville_massacre'),
  },
  {
    countryCode: 'ZAF', year: 1964, monthDay: '06-12',
    title: 'Rivonia Trial sentences Nelson Mandela to life imprisonment',
    summary:
      'Mandela and seven co-accused were convicted of sabotage and sentenced to life. His statement from the dock, ending with the words "an ideal for which I am prepared to die", was reported worldwide.',
    category: EventCategory.Politics, sourceUrl: wiki('Rivonia_Trial'),
  },
  {
    countryCode: 'ZAF', year: 1967, monthDay: '12-03',
    title: 'First human heart transplant',
    summary:
      'Christiaan Barnard performed the world\'s first human-to-human heart transplant at Groote Schuur Hospital in Cape Town. The patient, Louis Washkansky, survived 18 days before dying of pneumonia.',
    category: EventCategory.Science, sourceUrl: wiki('Christiaan_Barnard'),
  },
  {
    countryCode: 'ZAF', year: 1976, monthDay: '06-16',
    title: 'Soweto uprising',
    summary:
      'School students protesting the compulsory use of Afrikaans as a language of instruction were met with police gunfire, and unrest spread nationwide. The date is now commemorated as Youth Day.',
    category: EventCategory.Society, sourceUrl: wiki('Soweto_uprising'),
  },
  {
    countryCode: 'ZAF', year: 1990, monthDay: '02-11',
    title: 'Nelson Mandela released from prison',
    summary:
      'Mandela walked free after 27 years, nine days after President F. W. de Klerk unbanned the ANC and other organisations. Negotiations to end apartheid began that year.',
    category: EventCategory.Politics, sourceUrl: wiki('Nelson_Mandela'),
  },
  {
    countryCode: 'ZAF', year: 1994, monthDay: '04-27',
    title: 'First non-racial democratic election',
    summary:
      'South Africans of all races voted together for the first time, with queues stretching for kilometres. The ANC won and Nelson Mandela became president in May.',
    category: EventCategory.Politics, sourceUrl: wiki('1994_South_African_general_election'),
  },
  {
    countryCode: 'ZAF', year: 1996, monthDay: '12-10',
    title: 'Post-apartheid Constitution signed into law',
    summary:
      'The new constitution included an extensive Bill of Rights and created a Constitutional Court. It was the first national constitution in the world to prohibit discrimination on the basis of sexual orientation.',
    category: EventCategory.Politics, sourceUrl: wiki('Constitution_of_South_Africa'),
  },
  {
    countryCode: 'ZAF', year: 2010, monthDay: '06-11',
    title: 'FIFA World Cup opens in South Africa',
    summary:
      'South Africa became the first African country to host the tournament, playing across nine cities. Spain won the final in Johannesburg.',
    category: EventCategory.Culture, sourceUrl: wiki('2010_FIFA_World_Cup'),
  },

  // ----------------------------------------------------------------------- Egypt
  {
    countryCode: 'EGY', year: 1922, monthDay: '02-28',
    title: 'Britain recognises Egyptian independence',
    summary:
      'A unilateral declaration ended the protectorate and made Egypt a kingdom under Fuad I, though Britain retained control over defence, the Suez Canal and Sudan. Full sovereignty came only after 1952.',
    category: EventCategory.Politics, sourceUrl: wiki('Unilateral_Declaration_of_Egyptian_Independence'),
  },
  {
    countryCode: 'EGY', year: 1922, monthDay: '11-04',
    title: 'Tomb of Tutankhamun discovered',
    summary:
      'Howard Carter\'s team found the near-intact tomb of the boy pharaoh in the Valley of the Kings. The thousands of objects inside transformed Egyptology and set off a wave of popular "Egyptomania".',
    category: EventCategory.Culture, sourceUrl: wiki('Tutankhamun'),
  },
  {
    countryCode: 'EGY', year: 1952, monthDay: '07-23',
    title: 'Free Officers overthrow the monarchy',
    summary:
      'A coup led by the Free Officers deposed King Farouk, and a republic was declared the following year. Gamal Abdel Nasser emerged as leader and dominated Egyptian politics until 1970.',
    category: EventCategory.Politics, sourceUrl: wiki('Egyptian_revolution_of_1952'),
  },
  {
    countryCode: 'EGY', year: 1956, monthDay: '07-26',
    title: 'Nationalisation of the Suez Canal',
    summary:
      'Nasser nationalised the canal company to fund the Aswan High Dam after Western financing was withdrawn. Britain, France and Israel invaded in October but withdrew under US and Soviet pressure.',
    category: EventCategory.Politics, sourceUrl: wiki('Suez_Crisis'),
  },
  {
    countryCode: 'EGY', year: 1967, monthDay: '06-05',
    title: 'Six-Day War',
    summary:
      'Israel launched pre-emptive strikes that destroyed much of the Egyptian air force on the ground, and captured the Sinai Peninsula and Gaza Strip within six days. The defeat reshaped politics across the Arab world.',
    category: EventCategory.Conflict, sourceUrl: wiki('Six-Day_War'),
  },
  {
    countryCode: 'EGY', year: 1970, monthDay: '07-21',
    title: 'Aswan High Dam completed',
    summary:
      'The dam gave Egypt control of the Nile flood, added hydroelectric capacity and created Lake Nasser. It also displaced around 100,000 people and required the relocation of the Abu Simbel temples.',
    category: EventCategory.Science, sourceUrl: wiki('Aswan_Dam'),
  },
  {
    countryCode: 'EGY', year: 1973, monthDay: '10-06',
    title: 'October War begins',
    summary:
      'Egyptian forces crossed the Suez Canal and breached the Bar Lev Line on the Jewish holiday of Yom Kippur. Although the war ended inconclusively, it restored Egyptian confidence and opened the way to negotiations.',
    category: EventCategory.Conflict, sourceUrl: wiki('Yom_Kippur_War'),
  },
  {
    countryCode: 'EGY', year: 1978, monthDay: '09-17',
    title: 'Camp David Accords signed',
    summary:
      'Anwar Sadat and Menachem Begin signed frameworks for peace after thirteen days of talks hosted by Jimmy Carter. A formal Egypt–Israel treaty followed in 1979, and Sadat and Begin shared the Nobel Peace Prize.',
    category: EventCategory.Politics, sourceUrl: wiki('Camp_David_Accords'),
  },
  {
    countryCode: 'EGY', year: 2011, monthDay: '01-25',
    title: 'Egyptian revolution begins',
    summary:
      'Mass demonstrations centred on Tahrir Square in Cairo led to President Hosni Mubarak resigning on 11 February after nearly thirty years in power. It was part of the wider wave of Arab uprisings.',
    category: EventCategory.Society, sourceUrl: wiki('Egyptian_revolution_of_2011'),
  },

  // ------------------------------------------------------------------- Australia
  {
    countryCode: 'AUS', year: 1901, monthDay: '01-01',
    title: 'Federation of the Australian colonies',
    summary:
      'Six British colonies federated as the Commonwealth of Australia under a new constitution. The first federal parliament met in Melbourne in May.',
    category: EventCategory.Politics, sourceUrl: wiki('Federation_of_Australia'),
  },
  {
    countryCode: 'AUS', year: 1915, monthDay: '04-25',
    title: 'Landing at Gallipoli',
    summary:
      'Australian and New Zealand troops landed on the Gallipoli peninsula as part of a failed Allied campaign against the Ottoman Empire. The date is commemorated as Anzac Day.',
    category: EventCategory.Conflict, sourceUrl: wiki('Landing_at_Anzac_Cove'),
  },
  {
    countryCode: 'AUS', year: 1927, monthDay: '05-09',
    title: 'Federal Parliament moves to Canberra',
    summary:
      'Parliament opened in the purpose-built capital, chosen as a compromise between Sydney and Melbourne and laid out to a design by Walter Burley Griffin and Marion Mahony Griffin.',
    category: EventCategory.Politics, sourceUrl: wiki('Canberra'),
  },
  {
    countryCode: 'AUS', year: 1932, monthDay: '03-19',
    title: 'Sydney Harbour Bridge opens',
    summary:
      'The steel arch bridge opened after eight years of construction during the Depression, employing thousands. It remains the widest long-span bridge in the world.',
    category: EventCategory.Culture, sourceUrl: wiki('Sydney_Harbour_Bridge'),
  },
  {
    countryCode: 'AUS', year: 1942, monthDay: '02-19',
    title: 'Bombing of Darwin',
    summary:
      'Japanese aircraft attacked Darwin in two raids, killing at least 235 people in the largest attack ever mounted on Australia. It shifted Australian strategic reliance from Britain towards the United States.',
    category: EventCategory.Conflict, sourceUrl: wiki('Bombing_of_Darwin'),
  },
  {
    countryCode: 'AUS', year: 1967, monthDay: '05-27',
    title: 'Referendum on Aboriginal Australians passes',
    summary:
      'More than 90 percent voted to let the Commonwealth make laws for Aboriginal people and to include them in the census. It remains the largest "yes" vote in an Australian referendum.',
    category: EventCategory.Society, sourceUrl: wiki('1967_Australian_referendum_(Aboriginals)'),
  },
  {
    countryCode: 'AUS', year: 1973, monthDay: '10-20',
    title: 'Sydney Opera House opens',
    summary:
      'Queen Elizabeth II opened the building sixteen years after Jørn Utzon won the design competition, ten years late and far over budget. It became a UNESCO World Heritage Site in 2007.',
    category: EventCategory.Culture, sourceUrl: wiki('Sydney_Opera_House'),
  },
  {
    countryCode: 'AUS', year: 1992, monthDay: '06-03',
    title: 'Mabo decision recognises native title',
    summary:
      'The High Court held that the doctrine of terra nullius did not apply to Australia, recognising the land rights of the Meriam people. The Native Title Act followed in 1993.',
    category: EventCategory.Society, sourceUrl: wiki('Mabo_v_Queensland_(No_2)'),
  },
  {
    countryCode: 'AUS', year: 2000, monthDay: '09-15',
    title: 'Sydney Summer Olympics open',
    summary:
      'Australia hosted the Games for the second time, with Cathy Freeman lighting the cauldron and later winning the 400 metres. The organisation was widely praised.',
    category: EventCategory.Culture, sourceUrl: wiki('2000_Summer_Olympics'),
  },
  {
    countryCode: 'AUS', year: 2008, monthDay: '02-13',
    title: 'National Apology to the Stolen Generations',
    summary:
      'Prime Minister Kevin Rudd apologised in Parliament to Aboriginal and Torres Strait Islander children removed from their families under past government policies. The speech was broadcast live nationwide.',
    category: EventCategory.Society, sourceUrl: wiki('Apology_to_Australia%27s_Indigenous_peoples'),
  },
];
