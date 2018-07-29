'''provides two modes of accesing csv files: either by line or by column'''

class MyLine:
    '''used as line in the line mode:
      line[i] returns i'th entry on the line (from 0, first column is not included),
      line[string] returns entry in the column named string
'''

    def __init__(self,values,csvHandle):
        self.csv = csvHandle
        self.list = values     

    def __getitem__(self,i):
        if isinstance( i, int ):        
            return self.list[i]
        else:
            return self.list[ self.csv.keys[i] ]

    def __len__(self):
        return len(self.list)

class MyCsv:
    '''line mode: 
      data[string] will return line which has first item on the line 'string
      data[int] returns i'th row (from 0, header not included)
    # is understood as comment
    '''
    def __init__(self,fileName):
        self.keys = {}     
        self.index = {}
        self.lines = []

        lines = open(fileName).readlines()
   
        j = 0
        for j in range(len(lines)):
            if lines[j].strip()[0] == '#':
              continue
            elif len(self.keys) == 0:
                s = lines[j].split(',')
                for i in range(len(s)-1):
                    self.keys[s[i+1].strip()] = i
            else:
                s = lines[j].split(',')
                name = s[0]
                values = [item.strip() for item in s[1:]]

                self.index[name] = len(self.lines)
                self.lines.append( MyLine(values, self) ) 

    def __getitem__(self,i):
        if isinstance( i, int ):        
            return self.lines[i]
        else:
            return self.lines[ self.index[i] ]

    def __len__(self):
        return len(self.lines)


class MyCsvCol:
    '''column mode (colums are python lists): 
      data[string] will return column which has 'string' in the header
      data[int] returns i'th column (from 0, header not included)
    # is understood as comment line
    '''
    def __init__(self,fileName):
        self.keys = {}
        self.cols = []

        lines = open(fileName).readlines()
   
        j = 0
        for j in range(len(lines)):
            if lines[j].strip()[0] == '#':
              continue
            elif len(self.keys) == 0:
                s = lines[j].split(',')
                for i in range(len(s)):
                    self.keys[s[i].strip()] = i
                    self.cols.append( [] )
            else:
                s = lines[j].split(',')
                for i in range(len(s)):
                    item = s[i].strip()
                    try:
                      self.cols[i].append( item )
                    except IndexError:
                      print 'Error: line has length %d ' % (i+1)
                      print lines[j].strip()
                      raise IndexError

    def __getitem__(self,i):
        if isinstance( i, int ):        
            return self.cols[i]
        else:
            return self.cols[ self.keys[i] ]

    def __len__(self):
        return len(self.cols)






