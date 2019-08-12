from random import choice as choose
from random import random as r
from random import randint

CHANCE_TO_EXIT_AFTER_SENTENCE_TERMINATES = randint(25, 50)
DOUBLE_SPACE_BETWEEN_SENTENCES = True
MAXIMUM_TWEET_LENGTH = randint(randint(150, 275), 280)
PARAGRAPHS_TO_USE = randint(10, 30)
STRINGS_TO_DELETE_FROM_CORPUS = ["\n", "SUMMARY...", "DISCUSSION..."]
TWEETS_TO_GENERATE = 10


def chance(c):
    return c / 100. > r()


def terminatesSentence(word):
    if word[-1] in ["!", "?"]:
        return True
    # Period terminates, but not ellipses.
    elif word[-1] == ".":
        return word[-2] != "."
    else:
        return False


def addDictEntry(d, a, b):
    # If the key is not all uppercase (an abbreivation), make it lowercase.
    # Special exception for the article "a".
    # if (a == "A"):
    #     a = "a"
    # elif (a != None) and not (a.upper() == a):
    #     a = a.lower()

    try:
        d[a].append(b)
    except:
        d[a] = [b]


def superstrip(a):
    if a == None:
        return None
    simple = a
    for garbage in [".", "...", "?", "!", "(", ")", ","]:
        simple = simple.replace(garbage, "")
    return simple.lower()


def addAliasEntry(d, a):
    s = superstrip(a)
    try:
        if a not in d[s]:
            d[s].append(a)
    except:
        d[s] = [a]


paragraphs = []

with open("C:/Users/alpha/Documents/python/markov/spctext.txt") as f:
    thisParagraph = ""
    for line in f:
        for term in STRINGS_TO_DELETE_FROM_CORPUS:
            line = line.replace(term, "")

        if len(line.strip()) == 0:
            paragraphs.append(thisParagraph[:-1])
            thisParagraph = ""
        else:
            thisParagraph += line + " "
    if len(thisParagraph.strip()) > 0:
        paragraphs.append(thisParagraph[:-1])
    f.close()

d = dict()
# alias = dict()
# null = False
firstPara = randint(0, len(paragraphs) - PARAGRAPHS_TO_USE)

for para in paragraphs[firstPara:(firstPara + PARAGRAPHS_TO_USE)]:
    words = para.split(" ")
    while True:
        try:
            words.remove("")
        except:
            break

    # First word succeeds null.
    addDictEntry(d, None, words[0])

    for w in range(len(words) - 1):
        a, b = words[w:w + 2]
        if terminatesSentence(a):
            addDictEntry(d, a, None)
            addDictEntry(d, None, b)
        else:
            addDictEntry(d, a, b)

        # addAliasEntry(alias, a)

    # Last word preceeds null.
    addDictEntry(d, words[-1], None)
    # addAliasEntry(alias, words[-1])

# alias[None] = [None]

tweets = []
for t in range(TWEETS_TO_GENERATE):
    tweet = ""
    lastWord = None
    sentences = 0

    while len(tweet) < MAXIMUM_TWEET_LENGTH:
        # choices = []
        # for al in alias[superstrip(lastWord)]:
        #     choices += d[al]

        w = choose(d[lastWord])
        lastWord = w
        if w is None:
            if chance(CHANCE_TO_EXIT_AFTER_SENTENCE_TERMINATES):
                break
            tweet += " " if DOUBLE_SPACE_BETWEEN_SENTENCES else ""
            sentences += 1
        else:
            if len(tweet + w) > (MAXIMUM_TWEET_LENGTH - 3):
                tweet = tweet[:-1] + "..."
                break
            else:
                tweet += w + " "
    print(tweet)
    print()
    tweets.append(tweet)
